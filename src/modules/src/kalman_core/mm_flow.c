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

#include "mm_flow.h"
#include "log.h"
#include "platform_defaults.h"
#include "param.h"
#include <math.h>

#define FLOW_RESOLUTION 0.10f //We get the measurements in 10x the motion pixels (experimentally measured)

// Student-t robust weight on the flow update. The PMW3901 emits transient GARBAGE motion spikes
// when it blinds in the dark (flight 190829: deltaX/Y -> 8-26 px, driving vy to ~1 m/s); the
// light-aware std (flowdeck_v1v2.c) raises the nominal std, but a big transient spike can still
// leak through. This rejects it the gate-free way -- the same Student-t robust weight the wall /
// surface fusion uses: w = (nu+1)/(nu + innov^2/S), R inflated by 1/clip(w, wMin, 1). Normal flow
// (small innovation) -> w~1, untouched; a dark spike -> huge innovation -> w->wMin, rejected.
static uint8_t flowRobust = 1;   // 1 = Student-t robust weight on the flow update (default on)
static float flowNu = 5.0f;      // Student-t dof (flow_dof): heavy tails -> soft outlier rejection
static float flowWMin = 0.001f;  // robust weight floor: max R inflation = 1/wMin

// TODO remove the temporary test variables (used for logging)
static float predictedNX;
static float predictedNY;
static float measuredNX;
static float measuredNY;

static Axis3f flowdeckPos = { .axis = { FLOWDECK_POS_X, FLOWDECK_POS_Y, FLOWDECK_POS_Z } }; // In body coordinate system

// Inflation factor 1/sqrt(w) for the flow measurement std, from the Student-t robust weight.
// The flow H is sparse (only KC_STATE_Z and the velocity axis velIdx are non-zero), so the
// innovation variance S = H P H' + R is the 2-term quadratic below. Returns 1.0 when disabled or
// degenerate (no change to the stock update).
static float flowRobustStdScale(const kalmanCoreData_t* this, const float* h, int velIdx,
                                float std, float innov)
{
  if (!flowRobust) {
    return 1.0f;
  }
  float R = std * std;
  float S = h[KC_STATE_Z]  * h[KC_STATE_Z]  * this->P[KC_STATE_Z][KC_STATE_Z]
          + 2.0f * h[KC_STATE_Z] * h[velIdx] * this->P[KC_STATE_Z][velIdx]
          + h[velIdx] * h[velIdx] * this->P[velIdx][velIdx]
          + R;
  if (!(S > 0.0f)) {
    return 1.0f;
  }
  float w = (flowNu + 1.0f) / (flowNu + innov * innov / S);
  if (w > 1.0f) { w = 1.0f; }
  if (w < flowWMin) { w = flowWMin; }
  return 1.0f / sqrtf(w);  // R_eff = R/w  =>  std_eff = std / sqrt(w)
}

void kalmanCoreUpdateWithFlow(kalmanCoreData_t* this, const flowMeasurement_t *flow, const Axis3f *gyro)
{
  // Inclusion of flow measurements in the EKF done by two scalar updates

  // ~~~ Camera constants ~~~
  // The angle of aperture is guessed from the raw data register and thankfully look to be symmetric
  float Npix = 35.0;         // [pixels] (same in x and y)
  float thetapix = 0.71674f; // [rad] 2*sin(42/2); 42 degrees is the angle of aperture, here we computed the corresponding ground length

  //~~~ Extract states ~~~
  // Body rates
  float omegax_b = gyro->x * DEG_TO_RAD;
  float omegay_b = gyro->y * DEG_TO_RAD;
  float omegaz_b = gyro->z * DEG_TO_RAD;

  // Velocities in body frame
  float dx_b = this->S[KC_STATE_PX];
  float dy_b = this->S[KC_STATE_PY];

  // Height above origin
  float z_g = 0.0;

  // Saturate height in prediction and correction to avoid singularities
  if ( this->S[KC_STATE_Z] < 0.1f ) {
      z_g = 0.1;
  } else {
      z_g = this->S[KC_STATE_Z];
  }

  // Lever-arm induced translational velocity at camera
  // omega x r
  float v_cam_bx_add =  omegay_b * flowdeckPos.z - omegaz_b * flowdeckPos.y;
  float v_cam_by_add =  omegaz_b * flowdeckPos.x - omegax_b * flowdeckPos.z;
  
  // Effective camera point velocities in body frame
  float v_cam_bx = dx_b + v_cam_bx_add;
  float v_cam_by = dy_b + v_cam_by_add;

  // X velocity prediction and update
  // predicts the number of accumulated pixels in the x-direction
  float hx[KC_STATE_DIM] = {0};
  arm_matrix_instance_f32 Hx = {1, KC_STATE_DIM, hx};
  predictedNX = (flow->dt * Npix / thetapix) * ((v_cam_bx * this->R[2][2] / z_g) - omegay_b);
  measuredNX = flow->dpixelx*FLOW_RESOLUTION;

  // derive measurement equation with respect to dx (and z?)
  hx[KC_STATE_Z]  = (Npix * flow->dt / thetapix) * ((this->R[2][2] * v_cam_bx) / (-z_g * z_g));
  hx[KC_STATE_PX] = (Npix * flow->dt / thetapix) * (this->R[2][2] / z_g);

  //First update (robust: a dark garbage spike is rejected via the Student-t weight)
  float innovNX = measuredNX - predictedNX;
  float stdNX = flow->stdDevX * FLOW_RESOLUTION;
  kalmanCoreScalarUpdate(this, &Hx, innovNX, stdNX * flowRobustStdScale(this, hx, KC_STATE_PX, stdNX, innovNX));

  // Y velocity prediction and update
  float hy[KC_STATE_DIM] = {0};
  arm_matrix_instance_f32 Hy = {1, KC_STATE_DIM, hy};
  predictedNY = (flow->dt * Npix / thetapix ) * ((v_cam_by * this->R[2][2] / z_g) + omegax_b);
  measuredNY = flow->dpixely*FLOW_RESOLUTION;

  // derive measurement equation with respect to dy (and z?)
  hy[KC_STATE_Z]  = (Npix * flow->dt / thetapix) * ((this->R[2][2] * v_cam_by) / (-z_g * z_g));
  hy[KC_STATE_PY] = (Npix * flow->dt / thetapix) * (this->R[2][2] / z_g);

  // Second update (robust)
  float innovNY = measuredNY - predictedNY;
  float stdNY = flow->stdDevY * FLOW_RESOLUTION;
  kalmanCoreScalarUpdate(this, &Hy, innovNY, stdNY * flowRobustStdScale(this, hy, KC_STATE_PY, stdNY, innovNY));
}

/**
 * Predicted and measured values of the X and Y direction of the flowdeck
 */
LOG_GROUP_START(kalman_pred)

/**
 * @brief Flow sensor predicted dx  [pixels/frame]
 * 
 *  note: rename to kalmanMM.flowX?
 */
  LOG_ADD(LOG_FLOAT, predNX, &predictedNX)
/**
 * @brief Flow sensor predicted dy  [pixels/frame]
 * 
 *  note: rename to kalmanMM.flowY?
 */
  LOG_ADD(LOG_FLOAT, predNY, &predictedNY)
/**
 * @brief Flow sensor measured dx  [pixels/frame]
 * 
 *  note: This is the same as motion.deltaX, so perhaps remove this?
 */
  LOG_ADD(LOG_FLOAT, measNX, &measuredNX)
/**
 * @brief Flow sensor measured dy  [pixels/frame]
 * 
 *  note: This is the same as motion.deltaY, so perhaps remove this?
 */
  LOG_ADD(LOG_FLOAT, measNY, &measuredNY)
LOG_GROUP_STOP(kalman_pred)

/**
 * Flowdeck properties
 */
PARAM_GROUP_START(flowdeck)
  /**
   * @brief Flow deck position X (in meters, body frame)
   */
  PARAM_ADD_CORE(PARAM_FLOAT | PARAM_PERSISTENT, flowdeckPos_x, &flowdeckPos.x)
  /**
   * @brief Flow deck position Y (in meters, body frame)
   */
  PARAM_ADD_CORE(PARAM_FLOAT | PARAM_PERSISTENT, flowdeckPos_y, &flowdeckPos.y)
  /**
   * @brief Flow deck position Z (in meters, body frame)
   */
  PARAM_ADD_CORE(PARAM_FLOAT | PARAM_PERSISTENT, flowdeckPos_z, &flowdeckPos.z)
PARAM_GROUP_STOP(flowdeck)

/**
 * Flow update robustness: Student-t weight that rejects the PMW3901's dark garbage spikes.
 */
PARAM_GROUP_START(flowrob)
/**
 * @brief 1 = Student-t robust weight on the flow update (default), 0 = plain update.
 */
PARAM_ADD(PARAM_UINT8, enable, &flowRobust)
/**
 * @brief Student-t degrees of freedom: lower => softer rejection of a large flow innovation.
 */
PARAM_ADD(PARAM_FLOAT, nu, &flowNu)
/**
 * @brief Robust weight floor: the most the flow R can be inflated is 1/wMin.
 */
PARAM_ADD(PARAM_FLOAT, wMin, &flowWMin)
PARAM_GROUP_STOP(flowrob)
