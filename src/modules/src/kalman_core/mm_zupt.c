/**
 * ,---------,       ____  _ __
 * |  ,-^-,  |      / __ )(_) /_______________ _____  ___
 * | (  O  ) |     / __  / / __/ ___/ ___/ __ `/_  / / _ \
 * | / ,--'  |    / /_/ / / /_/ /__/ /  / /_/ / / /_/  __/
 *    +------`   /_____/_/\__/\___/_/   \__,_/ /___/\___/
 *
 * Crazyflie control firmware
 *
 * Copyright (C) 2021 Bitcraze AB
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
 */

#include "mm_zupt.h"
#include "param.h"

// Zero-velocity update (ZUPT) -- the principled replacement for the pre-climb estimator reset.
//
// While the craft sits on the ground it is stationary, so its body lateral velocity is exactly 0.
// The on-ground predict branch (kalman_core.c) integrates the raw accelerometer, and the 9-state
// EKF has NO accel-bias state, so the residual (accel bias + gravity that did not fully cancel
// because the attitude estimate is slightly off) is a roughly CONSTANT phantom horizontal
// acceleration -- double-integrated into a growing phantom velocity AND position, in a consistent
// direction. The estimator's only previous defence was to RESET (zero everything) at the last
// instant before the climb; any delay between that reset and takeoff re-opened the drift window
// and the controller flew the phantom offset into a diagonal drive.
//
// Fusing PX = PY = 0 each predict cycle while !quadIsFlying continuously cancels that drift, so the
// estimate stays clean on the ground for as long as the craft sits there -- no reset needed, and
// takeoff becomes timing-INDEPENDENT (a settle, a long preflight, etc. are all harmless). As a
// bonus it also holds X/Y at the origin through the early climb, so the wall references seed clean.
//
// Only PX/PY are pinned (the documented lateral-drift axes); Z is held by the down-ranger. The
// caller gates on !quadIsFlying so this can never fight real in-flight velocity, and even if
// supervisorIsFlying lags into the early climb, pinning lateral velocity to 0 there only helps the
// craft rise straight (there is no commanded lateral motion at takeoff).
static float zuptStdMs = 0.02f;  // [m/s] stationary-velocity std (tight: the craft is truly parked)
static uint8_t zuptEnable = 1;   // 1 = ZUPT on the ground (default on)

void kalmanCoreUpdateWithZupt(kalmanCoreData_t* this)
{
  if (!zuptEnable) {
    return;
  }
  float h[KC_STATE_DIM] = {0};
  arm_matrix_instance_f32 H = {1, KC_STATE_DIM, h};

  // measured body velocity = 0 on PX, then PY (sequential scalar updates; clear h between them).
  h[KC_STATE_PX] = 1.0f;
  kalmanCoreScalarUpdate(this, &H, -this->S[KC_STATE_PX], zuptStdMs);
  h[KC_STATE_PX] = 0.0f;

  h[KC_STATE_PY] = 1.0f;
  kalmanCoreScalarUpdate(this, &H, -this->S[KC_STATE_PY], zuptStdMs);
}

/**
 * Ground zero-velocity update (ZUPT) tuning -- the stationary-on-the-ground prior that replaces the
 * pre-climb estimator reset.
 */
PARAM_GROUP_START(zupt)
/**
 * @brief 1 = pin body lateral velocity to 0 while on the ground (default), 0 = off (revert to the
 * reset-only behaviour).
 */
PARAM_ADD(PARAM_UINT8, enable, &zuptEnable)
/**
 * @brief Stationary-velocity measurement std [m/s]. Tight, because a parked craft's velocity is
 * exactly 0; smaller = stronger pull against the accel-bias drift.
 */
PARAM_ADD(PARAM_FLOAT, std, &zuptStdMs)
PARAM_GROUP_STOP(zupt)
