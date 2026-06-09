/*
 *    ||          ____  _ __
 * +------+      / __ )(_) /_______________ _____  ___
 * | 0xBC |     / __  / / __/ ___/ ___/ __ `/_  / / _ \
 * +------+    / /_/ / / /_/ /__/ /  / /_/ / / /_/  __/
 *  ||  ||    /_____/_/\__/\___/_/   \__,_/ /___/\___/
 *
 * Copyright 2021, Bitcraze AB
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Foobar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 */
/* multiranger.c: Multiranger deck driver */
#include "deck.h"
#include "param.h"

#define DEBUG_MODULE "MR"

#include "system.h"
#include "debug.h"
#include "log.h"
#include "pca95x4.h"
#include "vl53l1x.h"
#include "range.h"
#include "estimator.h"
#include "static_mem.h"

#include "i2cdev.h"

#include "FreeRTOS.h"
#include "task.h"

#include <stdlib.h>

static bool isInit = false;
static bool isTested = false;
static bool isPassed = false;
static uint16_t filterMask = 1 << VL53L1_RANGESTATUS_RANGE_VALID;

// Fuse the UP beam into the estimator as the opposing-surface height anchor (mm_tof_surface).
// Default 0 = stock (UP beam log-only). Only valid returns within MR_CEILING_OUTLIER_MM are enqueued,
// so a missing deck / bad return simply feeds nothing (graceful fallback to down-beam-only).
static uint8_t mrUpFuse = 0;
static float mrUpStd = 0.05f;       // [m] fallback up-beam height std (used when no per-shot sigma)
#define MR_CEILING_OUTLIER_MM 5000  // reject >5 m / the 32767 censored sentinel

// Fuse the horizontal beams (front/back/left/right) into the estimator as the wall
// absolute-position anchor (mm_tof_walls). Default 0 = stock (log-only). Only valid returns within
// MR_CEILING_OUTLIER_MM are enqueued, so a wall out of range / a missing return simply feeds nothing
// (a room with only 2-3 walls in range still anchors each axis from whatever wall is present).
// Independent of mrUp.fuse (height) -- this owns horizontal X/Y.
static uint8_t mrWallsFuse = 0;
static float mrWallsStd = 0.02f;    // [m] fallback wall-range std (used when no per-shot sigma)

#define MR_PIN_UP PCA95X4_P0
#define MR_PIN_FRONT PCA95X4_P4
#define MR_PIN_BACK PCA95X4_P1
#define MR_PIN_LEFT PCA95X4_P6
#define MR_PIN_RIGHT PCA95X4_P2

NO_DMA_CCM_SAFE_ZERO_INIT static VL53L1_Dev_t devFront;
NO_DMA_CCM_SAFE_ZERO_INIT static VL53L1_Dev_t devBack;
NO_DMA_CCM_SAFE_ZERO_INIT static VL53L1_Dev_t devUp;
NO_DMA_CCM_SAFE_ZERO_INIT static VL53L1_Dev_t devLeft;
NO_DMA_CCM_SAFE_ZERO_INIT static VL53L1_Dev_t devRight;

// Per-beam VL53L1x return quality, captured from the ranging struct we already read
// every cycle (only RangeMilliMeter was kept before; the rest was discarded). These
// let the off-board fusion filter weight each beam by measured photon statistics
// (sigma, signal rate, ambient, status) instead of a geometric proxy, and recover the
// UNcensored range for grazing/oblique hits that range.* replaces with the 32767
// sentinel. Ranging config is untouched -- this only publishes discarded data.
typedef struct {
  float sigma;     // SigmaMilliMeter (mm): sensor's own range-std estimate
  float signal;    // SignalRateRtnMegaCps (MCPS): return strength / reflectance x geom
  float ambient;   // AmbientRateRtnMegaCps (MCPS): background light (for SNR)
  int16_t range;   // RangeMilliMeter (mm): UNcensored (range.* applies filterMask)
  uint8_t status;  // RangeStatus: 0=valid, else sigma/signal/wrap fail
} mrQuality_t;

static mrQuality_t qFront, qBack, qUp, qLeft, qRight;

static bool mrInitSensor(VL53L1_Dev_t *pdev, uint32_t pca95pin, char *name)
{
    bool status;

    // Bring up VL53 by releasing XSHUT
    pca95x4SetOutput(PCA95X4_DEFAULT_ADDRESS, pca95pin);
    // Let VL53 boot
    vTaskDelay(M2T(2));
    // Init VL53
    if (vl53l1xInit(pdev, I2C1_DEV))
    {
        DEBUG_PRINT("Init %s sensor [OK]\n", name);
        status = true;
    }
    else
    {
        DEBUG_PRINT("Init %s sensor [FAIL]\n", name);
        status = false;
    }

    return status;
}

static uint16_t mrGetMeasurementAndRestart(VL53L1_Dev_t *dev, mrQuality_t *q)
{
    VL53L1_Error status = VL53L1_ERROR_NONE;
    VL53L1_RangingMeasurementData_t rangingData;
    uint8_t dataReady = 0;
    uint16_t range;

    while (dataReady == 0)
    {
        status = VL53L1_GetMeasurementDataReady(dev, &dataReady);
        vTaskDelay(M2T(1));
    }

    status = VL53L1_GetRangingMeasurementData(dev, &rangingData);

    // Publish the return quality we would otherwise throw away (FixPoint16.16 -> float).
    // raw range is UNcensored so the filter can still use sigma-/signal-fail hits.
    q->sigma = rangingData.SigmaMilliMeter / 65536.0f;
    q->signal = rangingData.SignalRateRtnMegaCps / 65536.0f;
    q->ambient = rangingData.AmbientRateRtnMegaCps / 65536.0f;
    q->range = rangingData.RangeMilliMeter;
    q->status = rangingData.RangeStatus;

    if (filterMask & (1 << rangingData.RangeStatus))
    {
        range = rangingData.RangeMilliMeter;
    }
    else
    {
        range = 32767;
    }

    VL53L1_StopMeasurement(dev);
    status = VL53L1_StartMeasurement(dev);
    status = status;

    return range;
}

static void mrTask(void *param)
{
    VL53L1_Error status = VL53L1_ERROR_NONE;

    systemWaitStart();

    // Restart all sensors
    status = VL53L1_StopMeasurement(&devFront);
    status = VL53L1_StartMeasurement(&devFront);
    status = VL53L1_StopMeasurement(&devBack);
    status = VL53L1_StartMeasurement(&devBack);
    status = VL53L1_StopMeasurement(&devUp);
    status = VL53L1_StartMeasurement(&devUp);
    status = VL53L1_StopMeasurement(&devLeft);
    status = VL53L1_StartMeasurement(&devLeft);
    status = VL53L1_StopMeasurement(&devRight);
    status = VL53L1_StartMeasurement(&devRight);
    status = status;

    TickType_t lastWakeTime = xTaskGetTickCount();

    while (1)
    {
        vTaskDelayUntil(&lastWakeTime, M2T(100));
        uint16_t frontMm = mrGetMeasurementAndRestart(&devFront, &qFront);
        rangeSet(rangeFront, frontMm / 1000.0f);
        // Wall FRONT-beam anchor: feed a VALID front return to the estimator (mm_tof_walls -> X).
        // Same heteroscedastic pattern as the up beam: pass the per-shot VL53L1x sigma [m] as stdDev
        // so the handler trusts a clean return and self-weakens a far / grazing one. mrWallsStd is
        // the fallback when the sensor reports no sigma.
        if (mrWallsFuse != 0 && frontMm < MR_CEILING_OUTLIER_MM) {
          tofMeasurement_t front;
          front.timestamp = xTaskGetTickCount();
          front.distance = frontMm / 1000.0f;
          front.stdDev = (qFront.sigma > 0.0f) ? (qFront.sigma * 0.001f) : mrWallsStd;
          front.signal = qFront.signal;   // return strength (MCPS) + status for the world-plane quality gate
          front.status = qFront.status;   // (a grazing / edge / furniture hit self-reports here)
          estimatorEnqueueWallFront(&front);
        }
        uint16_t backMm = mrGetMeasurementAndRestart(&devBack, &qBack);
        rangeSet(rangeBack, backMm / 1000.0f);
        if (mrWallsFuse != 0 && backMm < MR_CEILING_OUTLIER_MM) {
          tofMeasurement_t back;
          back.timestamp = xTaskGetTickCount();
          back.distance = backMm / 1000.0f;
          back.stdDev = (qBack.sigma > 0.0f) ? (qBack.sigma * 0.001f) : mrWallsStd;
          back.signal = qBack.signal;
          back.status = qBack.status;
          estimatorEnqueueWallBack(&back);
        }
        uint16_t upMm = mrGetMeasurementAndRestart(&devUp, &qUp);
        rangeSet(rangeUp, upMm / 1000.0f);
        // Opposing-surface UP-beam anchor: feed a VALID up return to the estimator (mm_tof_surface).
        // Pass the per-shot VL53L1x sigma [m] as stdDev so the handler weights the up beam by its
        // OWN reported quality (heteroscedastic R) -- a far / weak return self-reports high sigma
        // and is trusted less. mrUpStd is the fallback when the sensor reports no sigma.
        if (mrUpFuse != 0 && upMm < MR_CEILING_OUTLIER_MM) {
          tofMeasurement_t up;
          up.timestamp = xTaskGetTickCount();
          up.distance = upMm / 1000.0f;
          up.stdDev = (qUp.sigma > 0.0f) ? (qUp.sigma * 0.001f) : mrUpStd; // qUp.sigma is [mm]
          estimatorEnqueueTOFSurfaceUp(&up);
        }
        uint16_t leftMm = mrGetMeasurementAndRestart(&devLeft, &qLeft);
        rangeSet(rangeLeft, leftMm / 1000.0f);
        // Wall LEFT-beam anchor: feed a VALID left return to the estimator (mm_tof_walls -> Y).
        if (mrWallsFuse != 0 && leftMm < MR_CEILING_OUTLIER_MM) {
          tofMeasurement_t left;
          left.timestamp = xTaskGetTickCount();
          left.distance = leftMm / 1000.0f;
          left.stdDev = (qLeft.sigma > 0.0f) ? (qLeft.sigma * 0.001f) : mrWallsStd;
          left.signal = qLeft.signal;
          left.status = qLeft.status;
          estimatorEnqueueWallLeft(&left);
        }
        uint16_t rightMm = mrGetMeasurementAndRestart(&devRight, &qRight);
        rangeSet(rangeRight, rightMm / 1000.0f);
        if (mrWallsFuse != 0 && rightMm < MR_CEILING_OUTLIER_MM) {
          tofMeasurement_t right;
          right.timestamp = xTaskGetTickCount();
          right.distance = rightMm / 1000.0f;
          right.stdDev = (qRight.sigma > 0.0f) ? (qRight.sigma * 0.001f) : mrWallsStd;
          right.signal = qRight.signal;
          right.status = qRight.status;
          estimatorEnqueueWallRight(&right);
        }
    }
}

static void mrInit()
{
    if (isInit)
    {
        return;
    }

    pca95x4Init();

    pca95x4ConfigOutput(PCA95X4_DEFAULT_ADDRESS,
                        ~(MR_PIN_UP |
                          MR_PIN_RIGHT |
                          MR_PIN_LEFT |
                          MR_PIN_FRONT |
                          MR_PIN_BACK));

    pca95x4ClearOutput(PCA95X4_DEFAULT_ADDRESS,
                       MR_PIN_UP |
                       MR_PIN_RIGHT |
                       MR_PIN_LEFT |
                       MR_PIN_FRONT |
                       MR_PIN_BACK);

    isInit = true;

    xTaskCreate(mrTask, MULTIRANGER_TASK_NAME, MULTIRANGER_TASK_STACKSIZE, NULL,
                MULTIRANGER_TASK_PRI, NULL);
}

static bool mrTest()
{
    if (isTested)
    {
        return isPassed;
    }

    isPassed = isInit;

    isPassed &= mrInitSensor(&devFront, MR_PIN_FRONT, "front");
    isPassed &= mrInitSensor(&devBack, MR_PIN_BACK, "back");
    isPassed &= mrInitSensor(&devUp, MR_PIN_UP, "up");
    isPassed &= mrInitSensor(&devLeft, MR_PIN_LEFT, "left");
    isPassed &= mrInitSensor(&devRight, MR_PIN_RIGHT, "right");

    isTested = true;

    return isPassed;
}

static const DeckDriver multiranger_deck = {
    .vid = 0xBC,
    .pid = 0x0C,
    .name = "bcMultiranger",

    .usedGpio = 0,
    .usedPeriph = DECK_USING_I2C,

    .init = mrInit,
    .test = mrTest,
};

DECK_DRIVER(multiranger_deck);

PARAM_GROUP_START(deck)

/**
 * @brief Nonzero if [Multi-ranger deck](%https://store.bitcraze.io/collections/decks/products/multi-ranger-deck) is attached
 */
PARAM_ADD_CORE(PARAM_UINT8 | PARAM_RONLY, bcMultiranger, &isInit)

PARAM_GROUP_STOP(deck)

PARAM_GROUP_START(multiranger)
/**
 * @brief Filter mask determining which range measurements is to be let through based on the range status of the VL53L1 chip
 */
PARAM_ADD(PARAM_UINT16, filterMask, &filterMask)

PARAM_GROUP_STOP(multiranger)

/**
 * Opposing-surface height: feed the UP beam to the estimator (mm_tof_surface).
 */
PARAM_GROUP_START(mrUp)
/**
 * @brief 1 = fuse the UP beam as the opposing-surface height anchor (mm_tof_surface), 0 = stock
 * (log-only). Pair with the down beam (zrange.surface=1) for furniture-robust height -- the up
 * surface holds height while a table occludes the down beam, and vice-versa for a lantern.
 */
PARAM_ADD(PARAM_UINT8, fuse, &mrUpFuse)
/**
 * @brief Fallback up-beam height measurement std [m] when the sensor reports no per-shot sigma.
 */
PARAM_ADD(PARAM_FLOAT, std, &mrUpStd)
PARAM_GROUP_STOP(mrUp)

/**
 * Wall absolute-position: feed the horizontal beams (front/back/left/right) to the estimator
 * (mm_tof_walls). front/back anchor world X, left/right anchor world Y.
 */
PARAM_GROUP_START(mrWalls)
/**
 * @brief 1 = fuse the four horizontal beams as the wall absolute-position anchor (mm_tof_walls),
 * 0 = stock (log-only). The wall references self-calibrate (a box pushed to a wall re-seats its
 * reference, X/Y hold), exactly as a table re-seats the down-surface height reference.
 */
PARAM_ADD(PARAM_UINT8, fuse, &mrWallsFuse)
/**
 * @brief Fallback wall-range measurement std [m] when the sensor reports no per-shot sigma.
 */
PARAM_ADD(PARAM_FLOAT, std, &mrWallsStd)
PARAM_GROUP_STOP(mrWalls)

/**
 * Per-beam VL53L1x return quality (front/back/left/right), sourced from the same
 * ranging struct as range.* but NOT censored by filterMask: 'raw*' carries the
 * grazing/oblique range that range.* drops to 32767. Lets the off-board fusion
 * filter weight each beam by measured photon statistics (sigma/signal/ambient/
 * status) instead of a geometric tilt/yaw proxy. Updated at the 10 Hz sensor rate.
 */
LOG_GROUP_START(mrq)
LOG_ADD(LOG_FLOAT, sigmaF, &qFront.sigma)
LOG_ADD(LOG_FLOAT, sigmaB, &qBack.sigma)
LOG_ADD(LOG_FLOAT, sigmaL, &qLeft.sigma)
LOG_ADD(LOG_FLOAT, sigmaR, &qRight.sigma)
LOG_ADD(LOG_FLOAT, signalF, &qFront.signal)
LOG_ADD(LOG_FLOAT, signalB, &qBack.signal)
LOG_ADD(LOG_FLOAT, signalL, &qLeft.signal)
LOG_ADD(LOG_FLOAT, signalR, &qRight.signal)
LOG_ADD(LOG_FLOAT, ambientF, &qFront.ambient)
LOG_ADD(LOG_FLOAT, ambientB, &qBack.ambient)
LOG_ADD(LOG_FLOAT, ambientL, &qLeft.ambient)
LOG_ADD(LOG_FLOAT, ambientR, &qRight.ambient)
LOG_ADD(LOG_INT16, rawF, &qFront.range)
LOG_ADD(LOG_INT16, rawB, &qBack.range)
LOG_ADD(LOG_INT16, rawL, &qLeft.range)
LOG_ADD(LOG_INT16, rawR, &qRight.range)
LOG_ADD(LOG_UINT8, statusF, &qFront.status)
LOG_ADD(LOG_UINT8, statusB, &qBack.status)
LOG_ADD(LOG_UINT8, statusL, &qLeft.status)
LOG_ADD(LOG_UINT8, statusR, &qRight.status)
// UP beam (ceiling) -- the opposing partner to the down/floor sensor for robust Z.
LOG_ADD(LOG_FLOAT, sigmaU, &qUp.sigma)
LOG_ADD(LOG_FLOAT, signalU, &qUp.signal)
LOG_ADD(LOG_INT16, rawU, &qUp.range)
LOG_ADD(LOG_UINT8, statusU, &qUp.status)
LOG_GROUP_STOP(mrq)
