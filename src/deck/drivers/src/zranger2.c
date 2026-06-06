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

// Runtime deweight of the onboard down-range fusion (analog of motion.flowStdFixed).
// 0 (default) => use the expStd model above (stock behaviour). >0 => OVERRIDE the enqueued
// std with this fixed value [m], so the onboard EKF barely trusts the raw down beam. This
// hands height authority to an external z (cf.extpos.send_extpos) for furniture-robust Z:
// a table sliding under steps the raw down beam to a tight ~2.5 mm-std reading that would
// otherwise weld stateEstimate.z to it (~400x our extpos trust) and drive the craft up.
static float downStdFixed = 0.0f;

// Range-RATE -> vertical-velocity fusion (furniture-immune fast damping). When velFuse != 0 the
// driver also enqueues d(range)/dt as a vertical velocity (mm_tof_rate). This keeps the fast
// vertical damping in the 1 kHz loop even when the ABSOLUTE down range is deweighted
// (downStdFixed) so an external z (extpos) owns height -- the fix for the deweight-causes-runaway
// failure. The rate is gated: a furniture EDGE is a one-tick jump (skipped); the surface gives
// the craft's true vertical speed.
static uint8_t surfaceFuse = 0;      // 1 -> route the down range through the opposing-surface
                                     //      height fusion (mm_tof_surface) instead of stock TOF
static uint8_t velFuse = 0;          // 1 -> fuse the range-rate as vertical velocity
static float velStd = 0.30f;         // [m/s] rate std ~ sqrt(2)*range_sigma/dt (1-sample diff noise)
static const float velDtMax = 0.10f; // [s] ignore stale gaps (re-seed instead of differencing)
static float ratePrevDist = 0.0f;    // [m] previous valid range
static uint32_t ratePrevTick = 0;
static bool rateHasPrev = false;

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
      if (downStdFixed > 0.0f) {
        stdDev = downStdFixed;  // runtime deweight: hand height authority to extpos-z
      }
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

      // Range-RATE -> vertical velocity (mm_tof_rate) for fast vertical damping. The furniture
      // defense is the opposing UP-beam surface (mm_tof_surface) + mm_tof_rate's own rate-vs-accel
      // sanity gate, so the driver no longer needs its own furniture step-gate here.
      if (velFuse != 0 && rateHasPrev) {
        float dt = (float)(now - ratePrevTick) * 0.001f; // 1 tick == 1 ms
        if (dt > 0.0f && dt < velDtMax) {
          tofMeasurement_t tofRate;
          tofRate.timestamp = now;
          tofRate.distance = (distance - ratePrevDist) / dt; // [m/s] packed into .distance
          tofRate.stdDev = velStd;
          estimatorEnqueueTOFRate(&tofRate);
        }
      }
      ratePrevDist = distance;
      ratePrevTick = now;
      rateHasPrev = true;
    } else {
      rateHasPrev = false; // outlier/gap -> re-seed the rate next valid sample (no diff across the gap)
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
 * @brief Fixed override [m] for the down-range measurement std fed to the EKF. 0 (default) =
 * use the built-in distance-dependent expStd model; >0 = deweight the onboard down beam so an
 * external z (cf.extpos.send_extpos) owns height (furniture-robust Z). Analog of motion.flowStdFixed.
 */
PARAM_ADD(PARAM_FLOAT, stdFixed, &downStdFixed)

/**
 * @brief 1 = also fuse the down-range RANGE-RATE as a vertical velocity (mm_tof_rate), 0 = off
 * (default, stock). Furniture-immune fast vertical damping: pair it with stdFixed>0 so the
 * absolute z is deweighted to an external extpos-z while the rate keeps the loop damped --
 * without it, deweighting alone runs the height away (no fast vertical sensor).
 */
PARAM_ADD(PARAM_UINT8, velFuse, &velFuse)
/**
 * @brief 1 = route the down range through the opposing-surface height fusion (mm_tof_surface:
 * gate-free Student-t, self-calibrating down reference that re-seats onto a table), 0 = stock TOF.
 * Pair with the up beam (mrUp.fuse=1) so the opposing surface holds height during an occlusion.
 * Keeps the fast onboard zranger as the height source.
 */
PARAM_ADD(PARAM_UINT8, surface, &surfaceFuse)
/**
 * @brief Measurement std [m/s] for the range-rate vertical-velocity fusion (velFuse).
 */
PARAM_ADD(PARAM_FLOAT, velStd, &velStd)

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
