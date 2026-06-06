/**
 *    ||          ____  _ __
 * +------+      / __ )(_) /_______________ _____  ___
 * | 0xBC |     / __  / / __/ ___/ ___/ __ `/_  / / _ \
 * +------+    / /_/ / / /_/ /__/ /  / /_/ / / /_/  __/
 *  ||  ||    /_____/_/\__/\___/_/   \__,_/ /___/\___/
 *
 * Crazyflie control firmware
 *
 * Copyright (C) 2021 BitCraze AB
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, in version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * vl53l0x.c: Time-of-flight distance sensor driver
 */

#define DEBUG_MODULE "ZR2"

#include "FreeRTOS.h"
#include "task.h"

#include "config.h"
#include "deck.h"
#include "system.h"
#include "debug.h"
#include "log.h"
#include "param.h"
#include "range.h"
#include "estimator.h"
#include "static_mem.h"

#include "i2cdev.h"
#include "zranger2.h"
#include "vl53l1x.h"

#include "cf_math.h"

// Measurement noise model
static const float expPointA = 2.5f;
static const float expStdA = 0.0025f; // STD at elevation expPointA [m]
static const float expPointB = 4.0f;
static const float expStdB = 0.2f;    // STD at elevation expPointB [m]
static float expCoeff;

// Furniture-robust HEIGHT is owned entirely by the opposing-surface fusion (mm_tof_surface):
// the floor (down) + ceiling (up) beams are fused as gate-free Student-t height with
// self-calibrating references. The down beam routes there when surfaceFuse != 0.
static uint8_t surfaceFuse = 0;      // 1 -> route the down range through the opposing-surface
                                     //      height fusion (mm_tof_surface) instead of stock TOF

#define RANGE_OUTLIER_LIMIT 5000 // the measured range is in [mm]

static uint16_t range_last = 0;

static bool isInit;

NO_DMA_CCM_SAFE_ZERO_INIT static VL53L1_Dev_t dev;

// FLOOR (down VL53L1x) return quality -- captured from the ranging struct already read
// each cycle and otherwise discarded. The floor signal/sigma is a direct read on the
// surface under the drone: a dark or featureless floor (which STARVES the PMW3901 flow
// deck) shows a weak/odd return here too, so this is an independent flow-health signal,
// and it doubles as the floor-vs-wall classifier for the horizontal beams. Ranging untouched.
static struct {
  float sigma;    // SigmaMilliMeter (mm)
  float signal;   // SignalRateRtnMegaCps (MCPS) -- floor reflectance/texture return
  float ambient;  // AmbientRateRtnMegaCps (MCPS)
  int16_t range;  // RangeMilliMeter (mm)
  uint8_t status; // RangeStatus
} zq;

static uint16_t zRanger2GetMeasurementAndRestart(VL53L1_Dev_t *dev)
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
    range = rangingData.RangeMilliMeter;

    // Publish the floor return quality we would otherwise throw away (FixPoint16.16 -> float).
    zq.sigma = rangingData.SigmaMilliMeter / 65536.0f;
    zq.signal = rangingData.SignalRateRtnMegaCps / 65536.0f;
    zq.ambient = rangingData.AmbientRateRtnMegaCps / 65536.0f;
    zq.range = rangingData.RangeMilliMeter;
    zq.status = rangingData.RangeStatus;

    VL53L1_StopMeasurement(dev);
    status = VL53L1_StartMeasurement(dev);
    status = status;

    return range;
}

void zRanger2Init(DeckInfo* info)
{
  if (isInit)
    return;

  if (vl53l1xInit(&dev, I2C1_DEV))
  {
      DEBUG_PRINT("Z-down sensor [OK]\n");
  }
  else
  {
    DEBUG_PRINT("Z-down sensor [FAIL]\n");
    return;
  }

  xTaskCreate(zRanger2Task, ZRANGER2_TASK_NAME, ZRANGER2_TASK_STACKSIZE, NULL, ZRANGER2_TASK_PRI, NULL);

  // pre-compute constant in the measurement noise model for kalman
  expCoeff = logf(expStdB / expStdA) / (expPointB - expPointA);

  isInit = true;
}

bool zRanger2Test(void)
{
  if (!isInit)
    return false;

  return true;
}

void zRanger2Task(void* arg)
{
  TickType_t lastWakeTime;

  systemWaitStart();

  // Restart sensor
  VL53L1_StopMeasurement(&dev);
  VL53L1_SetDistanceMode(&dev, VL53L1_DISTANCEMODE_MEDIUM);
  VL53L1_SetMeasurementTimingBudgetMicroSeconds(&dev, 25000);

  VL53L1_StartMeasurement(&dev);

  lastWakeTime = xTaskGetTickCount();

  while (1) {
    vTaskDelayUntil(&lastWakeTime, M2T(25));

    range_last = zRanger2GetMeasurementAndRestart(&dev);
    rangeSet(rangeDown, range_last / 1000.0f);

    // check if range is feasible and push into the estimator
    // the sensor should not be able to measure >5 [m], and outliers typically
    // occur as >8 [m] measurements
    if (range_last < RANGE_OUTLIER_LIMIT) {
      float distance = (float)range_last * 0.001f; // Scale from [mm] to [m]
      uint32_t now = xTaskGetTickCount();
      float stdDev = expStdA * (1.0f  + expf( expCoeff * (distance - expPointA)));
      if (surfaceFuse != 0) {
        // Opposing-surface DOWN-beam height fusion (mm_tof_surface): keeps the FAST onboard zranger
        // as the height source (no extpos-z latency), gate-free (Student-t), with a self-calibrating
        // down reference that re-seats onto a table. Pass the per-shot VL53L1x sigma [m] as stdDev
        // so the handler weights this beam by its OWN reported quality (heteroscedastic R); the
        // expStd model is the fallback when the sensor reports no sigma.
        tofMeasurement_t tofData;
        tofData.timestamp = now;
        tofData.distance = distance;
        tofData.stdDev = (zq.sigma > 0.0f) ? (zq.sigma * 0.001f) : stdDev; // zq.sigma is [mm]
        estimatorEnqueueTOFSurfaceDown(&tofData);
      } else {
        rangeEnqueueDownRangeInEstimator(distance, stdDev, now);
      }
    }
  }
}

static const DeckDriver zranger2_deck = {
  .vid = 0xBC,
  .pid = 0x0E,
  .name = "bcZRanger2",
  .usedGpio = 0,
  .usedPeriph = DECK_USING_I2C,

  .init = zRanger2Init,
  .test = zRanger2Test,
};

DECK_DRIVER(zranger2_deck);

PARAM_GROUP_START(deck)

/**
 * @brief Nonzero if [Z-ranger deck v2](%https://store.bitcraze.io/collections/decks/products/z-ranger-deck-v2) is attached
 */
PARAM_ADD_CORE(PARAM_UINT8 | PARAM_RONLY, bcZRanger2, &isInit)

PARAM_GROUP_STOP(deck)

/**
 * Down (floor) z-ranger fusion tuning.
 */
PARAM_GROUP_START(zrange)

/**
 * @brief 1 = route the down range through the opposing-surface height fusion (mm_tof_surface:
 * gate-free Student-t, self-calibrating down reference that re-seats onto a table), 0 = stock TOF.
 * Pair with the up beam (mrUp.fuse=1) so the opposing surface holds height during an occlusion.
 * Keeps the fast onboard zranger as the height source.
 */
PARAM_ADD(PARAM_UINT8, surface, &surfaceFuse)

PARAM_GROUP_STOP(zrange)

/**
 * Floor (down VL53L1x) return quality, sourced from the same ranging struct as
 * range.zrange. signal/sigma read the surface texture/reflectance under the drone -- an
 * independent flow-health signal (a flow-starving floor returns weak/odd here) and the
 * floor side of the floor-vs-wall classifier. Updated at the ~40 Hz zranger rate.
 */
LOG_GROUP_START(zrq)
LOG_ADD(LOG_FLOAT, sigma, &zq.sigma)
LOG_ADD(LOG_FLOAT, signal, &zq.signal)
LOG_ADD(LOG_FLOAT, ambient, &zq.ambient)
LOG_ADD(LOG_INT16, raw, &zq.range)
LOG_ADD(LOG_UINT8, status, &zq.status)
LOG_GROUP_STOP(zrq)
