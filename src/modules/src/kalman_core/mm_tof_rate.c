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

#include "mm_tof_rate.h"
#include "param.h"

// Sanity gate on |range-rate - accel-driven vz|: drop a rate that wildly disagrees with the
// accelerometer (a sensor glitch, or a furniture edge). GROUNDED at ~3-4*velStd (velStd~0.30 m/s
// is the 1-sample diff noise) so it never rejects a real climb/descent (where the accel moves PZ
// in step with the rate). The opposing CEILING absolute is the real furniture backstop now, so
// this is just a loose sanity bound, not the furniture defense it once had to be.
static float rateRobustGate = 1.0f; // [m/s]

// Fuse the down-ranger RANGE-RATE as a vertical-VELOCITY measurement on KC_STATE_PZ.
//
// Why: the absolute down distance (mm_tof) is what furniture under the craft fools -- a table
// welds the height estimate to a closer surface and drives the craft up. The RATE, however, is
// furniture-IMMUNE: a table edge is a one-tick spike the driver gates out, and on the table
// SURFACE the rate is still the craft's true vertical speed. Fusing the rate keeps the fast
// (sensor-rate) vertical DAMPING in the 1 kHz loop while an external absolute z
// (cf.extpos.send_extpos) owns the height -- the vertical analogue of the optical-flow deck,
// which lets the absolute zranger be deweighted (zrange.stdFixed) WITHOUT losing damping.
//
// Model (same geometry as mm_tof): d = z / cos(angle)  =>  d_dot = vz_world / cos(angle),
// with vz_world = (R * v_body)_z = R[2][0]*vx + R[2][1]*vy + R[2][2]*vz (body-frame velocities).
void kalmanCoreUpdateWithTofRate(kalmanCoreData_t* this, tofMeasurement_t *tofRate)
{
  float h[KC_STATE_DIM] = {0};
  arm_matrix_instance_f32 H = {1, KC_STATE_DIM, h};

  // Skip when near-horizontal (the model diverges as R[2][2] -> 0), exactly as mm_tof.
  if (fabsf(this->R[2][2]) > 0.1f && this->R[2][2] > 0.0f) {
    float angle = fabsf(acosf(this->R[2][2])) - DEG_TO_RAD * (15.0f / 2.0f);
    if (angle < 0.0f) {
      angle = 0.0f;
    }
    float cosAngle = cosf(angle);

    float predictedRate = (this->R[2][0] * this->S[KC_STATE_PX]
                         + this->R[2][1] * this->S[KC_STATE_PY]
                         + this->R[2][2] * this->S[KC_STATE_PZ]) / cosAngle;
    float measuredRate = tofRate->distance; // [m/s] -- driver packs the range rate into .distance
    float innov = measuredRate - predictedRate;

    // Robust gate: reject a rate that disagrees with the accel-driven vertical velocity -- a
    // gradual table slide (the per-tick range step is small enough to slip the absolute gate)
    // reads a descent the craft is NOT making, so vz must stay on the accelerometer.
    if (fabsf(innov) <= rateRobustGate) {
      h[KC_STATE_PX] = this->R[2][0] / cosAngle;
      h[KC_STATE_PY] = this->R[2][1] / cosAngle;
      h[KC_STATE_PZ] = this->R[2][2] / cosAngle;

      kalmanCoreScalarUpdate(this, &H, innov, tofRate->stdDev);
    }
  }
}

PARAM_GROUP_START(zrate)
/**
 * @brief |range-rate - accel-driven vz| [m/s] above which the rate is rejected as furniture (a
 * gradual table slide). The accelerometer is the anchor: it knows the craft is not descending.
 */
PARAM_ADD(PARAM_FLOAT, robustGate, &rateRobustGate)
PARAM_GROUP_STOP(zrate)
