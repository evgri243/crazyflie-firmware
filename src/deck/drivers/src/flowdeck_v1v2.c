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
/* flowdeck_v1v2.c: Flow deck driver */
#include "FreeRTOS.h"
#include "task.h"

#include "deck.h"
#include "debug.h"
#include "system.h"
#include "log.h"
#include "param.h"
#include "pmw3901.h"
#include "sleepus.h"

#include "stabilizer_types.h"
#include "estimator.h"
#include "estimator.h"

#include "cf_math.h"

#include "usec_time.h"
#include <stdlib.h>

#define AVERAGE_HISTORY_LENGTH 4
#define OULIER_LIMIT 100
#define LP_CONSTANT 0.8f
// #define USE_LP_FILTER
// #define USE_MA_SMOOTHING

#if defined(USE_MA_SMOOTHING)
static struct {
  float32_t averageX[AVERAGE_HISTORY_LENGTH];
  float32_t averageY[AVERAGE_HISTORY_LENGTH];
  size_t ptr;
} pixelAverages;
#endif

float dpixelx_previous = 0;
float dpixely_previous = 0;

static uint8_t outlierCount = 0;
static float stdFlow = 2.0f;

static bool isInit1 = false;
static bool isInit2 = false;

motionBurst_t currentMotion;

// Disables pushing the flow measurement in the EKF
static bool useFlowDisabled = false;

// Turn on adaptive standard deviation for the kalman filter
static bool useAdaptiveStd = false;

// Set standard deviation flow
// (will not work if useAdaptiveStd is on)
static float flowStdFixed = 2.0f;

// --- Light-aware flow std (honest R_flow) ---
// The PMW3901 BLINDS in the dark: it stops reporting ~0 and starts emitting GARBAGE motion spikes
// (flight 190829: deltaX/Y jumped to 8-26 px at squal 11 / contrast 7, driving vy to ~1 m/s). The
// fixed std fused that garbage straight into velocity. So ramp the flow std from base (good light)
// toward stdBlind (effectively ignore flow) as the sensor's OWN quality collapses -- keyed on squal
// (lock confidence) AND contrast = maxRawData - minRawData (texture/light), the cleaner blindness
// signals than shutter (which saturates at 8191). Anchors re-calibrated against real hover logs (see
// the static block below): trustworthy at squal>=20 / contrast>=60, blind at squal<=12 / contrast<=15.
// Default on; this is the flow's honest variance, so the EKF self-down-weights blind flow only when
// the sensor is truly near the garbage band, and the walls + accel carry velocity there.
static uint8_t useQualityStd = 1;        // 1 -> squal+contrast ramp (default); else adaptive/fixed
// Calibrated against real hover logs: in NORMAL flight this deck reads squal 18-42 (median ~28) --
// that is healthy flow, NOT blindness. The dark GARBAGE case (flight 190829) was squal 11-14. The
// first ramp (good=60/floor=15) put normal flight mid-ramp and downweighted ~7x for the WHOLE hover,
// softening velocity damping enough to let the craft drift (then the legacy wall refs follow it).
// So the ramp now trusts the normal range fully and only fades as squal approaches the garbage band.
// The Student-t robust weight (flowrob, mm_flow.c) stays the real catcher for transient dark spikes,
// so a gentle std ramp loses no dark protection.
static float flowSqualGood = 20.0f;      // squal at/above which flow is fully trusted (normal flight 18-42)
static float flowSqualFloor = 12.0f;     // squal at/below which flow is blind (garbage band ~11-14)
static float flowContrastGood = 60.0f;   // (maxRaw-minRaw) at/above which flow is fully trusted (sweep: dim~50, bright~180)
static float flowContrastFloor = 15.0f;  // (maxRaw-minRaw) at/below which flow is blind (dark~7)
static float flowStdBlind = 20.0f;       // flow std when fully blind (effective ~2 m/s -> ignored)

#define NCS_PIN DECK_GPIO_IO3


static void flowdeckTask(void *param)
{
  systemWaitStart();

  uint64_t lastTime  = usecTimestamp();
  while(1) {
    vTaskDelay(10);

    pmw3901ReadMotion(NCS_PIN, &currentMotion);

    // Flip motion information to comply with sensor mounting
    // (might need to be changed if mounted differently)
    int16_t accpx = -currentMotion.deltaY;
    int16_t accpy = -currentMotion.deltaX;

    // Outlier removal
    if (abs(accpx) < OULIER_LIMIT && abs(accpy) < OULIER_LIMIT) {

    if (useQualityStd)
    {
      // Honest R_flow: ramp the std from base (flowStdFixed) toward flowStdBlind as the sensor
      // quality collapses, so blind flow (which spikes garbage in the dark) self-down-weights to
      // nothing. squalScore/contrastScore each ramp 1->0 over good->floor; the WORST one dominates
      // (min) -- either signal saying "blind" suffices. q=1 -> base, q=0 -> blind.
      float squalScore = ((float)currentMotion.squal - flowSqualFloor) / (flowSqualGood - flowSqualFloor);
      float contrast = (float)currentMotion.maxRawData - (float)currentMotion.minRawData;
      float contrastScore = (contrast - flowContrastFloor) / (flowContrastGood - flowContrastFloor);
      if (squalScore < 0.0f) { squalScore = 0.0f; }
      if (squalScore > 1.0f) { squalScore = 1.0f; }
      if (contrastScore < 0.0f) { contrastScore = 0.0f; }
      if (contrastScore > 1.0f) { contrastScore = 1.0f; }
      float q = (squalScore < contrastScore) ? squalScore : contrastScore;
      stdFlow = flowStdFixed + (flowStdBlind - flowStdFixed) * (1.0f - q);
    }
    else if (useAdaptiveStd)
    {
      // The standard deviation is fitted by measurements flying over low and high texture
      //   and looking at the shutter time
      float shutter_f = (float)currentMotion.shutter;
      stdFlow=0.0007984f *shutter_f + 0.4335f;
      if (stdFlow < 0.1f) stdFlow=0.1f;
    } else {
      stdFlow = flowStdFixed;
    }


    // Form flow measurement struct and push into the EKF
    flowMeasurement_t flowData;
    flowData.stdDevX = stdFlow;
    flowData.stdDevY = stdFlow;
    flowData.dt = (float)(usecTimestamp()-lastTime)/1000000.0f;
    // we do want to update dt every measurement and not only in the ones with detected motion,
    // as we work with instantaneous gyro and velocity values in the update function
    // (meaning assuming the current measurements over all of dt)
    lastTime = usecTimestamp();



#if defined(USE_MA_SMOOTHING)
      // Use MA Smoothing
      pixelAverages.averageX[pixelAverages.ptr] = (float32_t)accpx;
      pixelAverages.averageY[pixelAverages.ptr] = (float32_t)accpy;

      float32_t meanX;
      float32_t meanY;

      arm_mean_f32(pixelAverages.averageX, AVERAGE_HISTORY_LENGTH, &meanX);
      arm_mean_f32(pixelAverages.averageY, AVERAGE_HISTORY_LENGTH, &meanY);

      pixelAverages.ptr = (pixelAverages.ptr + 1) % AVERAGE_HISTORY_LENGTH;

      flowData.dpixelx = (float)meanX;   // [pixels]
      flowData.dpixely = (float)meanY;   // [pixels]
#elif defined(USE_LP_FILTER)
      // Use LP filter measurements
      flowData.dpixelx = LP_CONSTANT * dpixelx_previous + (1.0f - LP_CONSTANT) * (float)accpx;
      flowData.dpixely = LP_CONSTANT * dpixely_previous + (1.0f - LP_CONSTANT) * (float)accpy;
      dpixelx_previous = flowData.dpixelx;
      dpixely_previous = flowData.dpixely;
#else
      // Use raw measurements
      flowData.dpixelx = (float)accpx;
      flowData.dpixely = (float)accpy;
#endif
      // Push measurements into the estimator if flow is not disabled
      //    and the PMW flow sensor indicates motion detection
      if (!useFlowDisabled && currentMotion.motion == 0xB0) {
        estimatorEnqueueFlow(&flowData);
      }
    } else {
      outlierCount++;
    }
  }
}

static void flowdeck1Init()
{
  if (isInit1 || isInit2) {
    return;
  }

  // Initialize the VL53L0 sensor using the zRanger deck driver
  const DeckDriver *zRanger = deckFindDriverByName("bcZRanger");
  zRanger->init(NULL);

  if (pmw3901Init(NCS_PIN))
  {
    xTaskCreate(flowdeckTask, FLOW_TASK_NAME, FLOW_TASK_STACKSIZE, NULL,
                FLOW_TASK_PRI, NULL);

    isInit1 = true;
  }
}

static bool flowdeck1Test()
{
  if (!isInit1) {
    DEBUG_PRINT("Error while initializing the PMW3901 sensor\n");
    return false;
  }

  // Test the VL53L0 driver
  const DeckDriver *zRanger = deckFindDriverByName("bcZRanger");

  return zRanger->test();
}

static const DeckDriver flowdeck1_deck = {
  .vid = 0xBC,
  .pid = 0x0A,
  .name = "bcFlow",
  .usedGpio = DECK_USING_IO_3,
  .usedPeriph = DECK_USING_I2C | DECK_USING_SPI,
  .requiredEstimator = StateEstimatorTypeKalman,

  .init = flowdeck1Init,
  .test = flowdeck1Test,
};

DECK_DRIVER(flowdeck1_deck);

static void flowdeck2Init()
{
  if (isInit1 || isInit2) {
    return;
  }

  // Initialize the VL53L1 sensor using the zRanger deck driver
  const DeckDriver *zRanger = deckFindDriverByName("bcZRanger2");
  zRanger->init(NULL);

  if (pmw3901Init(NCS_PIN))
  {
    xTaskCreate(flowdeckTask, FLOW_TASK_NAME, FLOW_TASK_STACKSIZE, NULL,
                FLOW_TASK_PRI, NULL);

    isInit2 = true;
  }
}

static bool flowdeck2Test()
{
  if (!isInit2) {
    DEBUG_PRINT("Error while initializing the PMW3901 sensor\n");
    return false;
  }

  // Test the VL53L1 driver
  const DeckDriver *zRanger = deckFindDriverByName("bcZRanger2");

  return zRanger->test();
}

static const DeckDriver flowdeck2_deck = {
  .vid = 0xBC,
  .pid = 0x0F,
  .name = "bcFlow2",

  .usedGpio = DECK_USING_IO_3,
  .usedPeriph = DECK_USING_I2C | DECK_USING_SPI,
  .requiredEstimator = StateEstimatorTypeKalman,

  .init = flowdeck2Init,
  .test = flowdeck2Test,
};

DECK_DRIVER(flowdeck2_deck);

/**
 * Logging variables of the motion sensor of the flowdeck
 */
LOG_GROUP_START(motion)
/**
 * @brief True if motion occurred since the last measurement
 */
LOG_ADD(LOG_UINT8, motion, &currentMotion.motion)
/**
 * @brief Flow X  measurement  [flow/fr]
 */
LOG_ADD(LOG_INT16, deltaX, &currentMotion.deltaX)
/**
 * @brief Flow Y measurement [flow/fr]
 */
LOG_ADD(LOG_INT16, deltaY, &currentMotion.deltaY)
/**
 * @brief Shutter time [clock cycles]
 */
LOG_ADD(LOG_UINT16, shutter, &currentMotion.shutter)
/**
 * @brief Maximum raw data value in frame
 */
LOG_ADD(LOG_UINT8, maxRaw, &currentMotion.maxRawData)
/**
 * @brief Minimum raw data value in frame
 */
LOG_ADD(LOG_UINT8, minRaw, &currentMotion.minRawData)
/**
 * @brief Average raw data value
 */
LOG_ADD(LOG_UINT8, Rawsum, &currentMotion.rawDataSum)
/**
 * @brief Counted flow outliers excluded from the estimator
 */
LOG_ADD(LOG_UINT8, outlierCount, &outlierCount)
/**
 * @brief Count of surface feature
 */
LOG_ADD(LOG_UINT8, squal, &currentMotion.squal)
/**
 * @brief Standard deviation of flow measurement
 */
LOG_ADD(LOG_FLOAT, std, &stdFlow)
LOG_GROUP_STOP(motion)

/**
 * Settings and parameters for handling of the flowdecks
 * measurements
 */
PARAM_GROUP_START(motion)
/**
 * @brief Nonzero to not push the flow measurement in the EKF (default: 0)
 */
PARAM_ADD(PARAM_UINT8, disable, &useFlowDisabled)
/**
 * @brief Nonzero to turn on adaptive standard deviation estimation (default: 0)
 */
PARAM_ADD(PARAM_UINT8, adaptive, &useAdaptiveStd)
/**
 * @brief Set standard deviation flow measurement (default: 2.0f)
 */
PARAM_ADD_CORE(PARAM_FLOAT | PARAM_PERSISTENT, flowStdFixed, &flowStdFixed)
/**
 * @brief 1 = light-aware flow std: ramp base->blind as squal/contrast collapse (honest R_flow,
 * default). Takes precedence over `adaptive`. 0 = use adaptive(shutter)/fixed.
 */
PARAM_ADD(PARAM_UINT8, qualityStd, &useQualityStd)
PARAM_ADD(PARAM_FLOAT, squalGood, &flowSqualGood)
PARAM_ADD(PARAM_FLOAT, squalFloor, &flowSqualFloor)
PARAM_ADD(PARAM_FLOAT, contrastGood, &flowContrastGood)
PARAM_ADD(PARAM_FLOAT, contrastFloor, &flowContrastFloor)
PARAM_ADD(PARAM_FLOAT, stdBlind, &flowStdBlind)
PARAM_GROUP_STOP(motion)

PARAM_GROUP_START(deck)

/**
 * @brief Nonzero if Flow deck v1 is attached
 */
PARAM_ADD_CORE(PARAM_UINT8 | PARAM_RONLY, bcFlow, &isInit1)

/**
 * @brief Nonzero if [Flow deck v2](%https://store.bitcraze.io/collections/decks/products/flow-deck-v2) is attached
 */
PARAM_ADD_CORE(PARAM_UINT8 | PARAM_RONLY, bcFlow2, &isInit2)

PARAM_GROUP_STOP(deck)
