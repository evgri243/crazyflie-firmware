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

#include "mm_tof_walls.h"
#include "param.h"
#include "log.h"
#include <math.h>

// ============================================================================
//  Wall absolute-position fusion (FRONT + BACK + LEFT + RIGHT horizontal beams)
//
//  The 4-wall horizontal mirror of mm_tof_surface.c. Where the vertical module
//  anchors HEIGHT against the floor (down beam) and ceiling (up beam), this module
//  anchors world X against the front/back beams and world Y against the left/right
//  beams. Each beam is a direct world-POSITION pseudo-measurement with a per-wall
//  self-calibrating reference, so a box pushed to a wall is just as step-like as a
//  table sliding under the down beam -- the reference absorbs it and X/Y do not move.
//
//  Naming: the references are named by the BEAM that sees them (front/back/left/right),
//  NOT by a fixed wall, because each reference SELF-CALIBRATES to whatever surface is
//  actually in the beam (the room wall, or a box pushed in front of it). The beam
//  direction is the only thing always true.
//
//  The whole machinery -- Student-t robust scalar update, jump-diffusion re-seat,
//  heteroscedastic R from the per-shot VL53L1x sigma, on-ground capture, the auxiliary
//  reference KF -- is IDENTICAL to mm_tof_surface.c. Only the geometry differs:
//
//   GEOMETRY (from the tested fusion/state.py SIGN table; the firmware uses the
//   rotation matrix R directly rather than roll/pitch angles).
//     beam   pos_axis      sign   beam world dir   dz (world-Z of beam)   cosSlant
//     front  KC_STATE_X    -1     +body X          dz =  R[2][0]          sqrt(1 - R[2][0]^2)
//     back   KC_STATE_X    +1     -body X          dz = -R[2][0]          sqrt(1 - R[2][0]^2)
//     left   KC_STATE_Y    -1     +body Y          dz =  R[2][1]          sqrt(1 - R[2][1]^2)
//     right  KC_STATE_Y    +1     -body Y          dz = -R[2][1]          sqrt(1 - R[2][1]^2)
//   The horizontal distance to the wall is r_h = range * cosSlant; the pseudo-measurement
//   is measured_pos = wallRef + sign * r_h, h[pos_axis] = 1 (the mm_tof_surface height
//   form with pos_axis in place of KC_STATE_Z).
//
//   TWO X/Y-specific defenses on top of the vertical module:
//   - RAY-CAST wall-vs-floor down-weight (fusion/sensors.py rangeabs_update): a beam that
//     descends (dz < 0) reaches the floor (world Z=0) at slant s_floor = -z/dz. If that is
//     nearer than the predicted wall slant, the beam is probably hitting the floor, not the
//     wall, so R is softly inflated by 1/w_wall (a sigmoid). dz >= 0 (level/up) => no floor.
//   - YAW-HOLD down-weight: the beam->axis mapping assumes the takeoff heading. Capture yaw0
//     while on the ground, then inflate R by 1 + (max(|roll|,|pitch|,|yaw_dev|)/tiltRScaleDeg)^2
//     so the walls fade out smoothly as the craft rotates (out of envelope by design).
// ============================================================================

// --- ported from the validated Python fusion config (fusion/state.py defaults) ---
//   These mirror the zsurf defaults 1:1 (the wall fusion shares the surface tuning);
//   the per-wall random-walk std is a touch faster than the stiff floor datum.
static float wallNu          = 4.0f;    // Student-t dof (abs_dof): heavy tails -> soft rejection
static float wallWMin        = 0.001f;  // robust weight floor (abs_w_min): max R inflation = 1/wMin
static float wallSigFloorM   = 0.005f;  // meas_sigma_floor_mm: VL53L1x sigma is floored here [m]
static float wallSigModelM   = 0.020f;  // meas_sigma_model_mm: wall non-planarity / mount [m]
static float wallStdFallbackM = 0.020f; // abs_meas_std_m: used when the per-shot sigma is invalid
static float wallJumpCouple  = 0.5f;    // abs_jump_couple: fraction of jumpSize^2 injected at w->0
static float wallJumpSizeM   = 0.35f;   // abs_jump_size_m: wall-swap step scale [m]
static float wallAVarMaxM2   = 0.25f;   // abs_a_var_max_m2: cap on a reference's variance [m^2]
static float wallReseatWMax  = 0.40f;   // reseat_w_max: only re-seat below this w (a clear swap)
static float wallRwStd       = 0.01f;   // per-wall reference random-walk std [m/sqrt(s)]
// ^ Honest static-wall prior. 0.10 asserted the wall wanders ~10 cm/s, which let the reference CHASE
//   the craft on a smooth translation (Kref ~0.83/step) -- the residual diagonal drift. A wall does
//   not move; its genuine slow change (sub-degree yaw bias, mount creep) is ~mm/s, so 0.01 (Kref
//   ~0.18). Sensor noise still lives in R (per-shot sigma + sigModel + tilt/floor inflation); a box
//   still re-seats via the Student-t jump (var injected up to aVarMax), independent of rw.
static float wallInitVarM2   = 0.01f;   // init_ref_std_m^2: reference variance at (re)seed [m^2]
static float wallTiltRScaleDeg = 20.0f; // tilt_r_scale_deg: yaw/tilt R-inflation scale [deg]
// ray-cast wall-vs-floor model (fusion/sensors.py rangeabs_update)
static float wallRaycastFloorScaleM = 0.15f; // abs_raycast_floor_scale_m: slant-margin sigmoid scale
static float wallRaycastCosMin      = 0.30f; // abs_raycast_cos_min: clamp on the horizontal fraction
static float wallRaycastFloorWMin   = 0.01f; // abs_raycast_floor_w_min: floor confidence lower clamp
static const float wallDtMaxS = 0.20f;  // clamp the random-walk dt across a sensor gap [s]

typedef struct {
  float ref;            // wall position along its axis [m] (~ preflight wall range at takeoff)
  float var;            // reference variance [m^2]
  uint32_t lastMs;      // tick of the last update (random-walk dt)
  bool seeded;          // reference initialised
  int posAxis;          // KC_STATE_X (front/back) or KC_STATE_Y (left/right)
  int rRow;             // rotation-matrix column index for dz/cosSlant: 0 for X, 1 for Y
  float sign;           // +1 (back/right: pos = ref - r_h) or -1 (front/left: pos = ref + r_h)
  float dzSign;         // multiplies R[2][rRow] to get the beam's world-Z (+1 front/left, -1 else)
  // diagnostics (first flights)
  float logRef;         // current reference [m]
  float logInnov;       // position innovation [m]
  float logW;           // robust weight (<1 => down-weighted; tiny => wall swap)
} wallRef_t;

// Module-static wall references (no Kalman states, exactly like the down/up surfaces). These persist
// across flights, so kalmanCoreWallReset() MUST be called on every estimator reset to clear `seeded`
// (otherwise a second flight inherits the previous flight's re-seated box/hand reference and flies on
// it -- flight 150720 drove into the table because fRef was stale at 0.35). Unlike the down beam,
// the horizontal beams do not reliably re-capture on the ground, so seeding happens at the first
// in-flight frame (X,Y ~ 0 at takeoff => ref ~ the true wall range); the reset makes every flight
// behave like the first post-flash flight (150319), which seeded correctly and held position.
//
// SIGN/dzSign per the state.py table: front pos = ref - r_h (sign -1, dz = +R[2][0]); back pos =
// ref + r_h (sign +1, dz = -R[2][0]); left pos = ref - r_h (sign -1, dz = +R[2][1]); right pos =
// ref + r_h (sign +1, dz = -R[2][1]).
static wallRef_t frontWall = {
  .ref = 0.0f, .var = 0.01f, .lastMs = 0, .seeded = false,
  .posAxis = KC_STATE_X, .rRow = 0, .sign = -1.0f, .dzSign = 1.0f,
  .logRef = 0.0f, .logInnov = 0.0f, .logW = 1.0f,
};
static wallRef_t backWall = {
  .ref = 0.0f, .var = 0.01f, .lastMs = 0, .seeded = false,
  .posAxis = KC_STATE_X, .rRow = 0, .sign = 1.0f, .dzSign = -1.0f,
  .logRef = 0.0f, .logInnov = 0.0f, .logW = 1.0f,
};
static wallRef_t leftWall = {
  .ref = 0.0f, .var = 0.01f, .lastMs = 0, .seeded = false,
  .posAxis = KC_STATE_Y, .rRow = 1, .sign = -1.0f, .dzSign = 1.0f,
  .logRef = 0.0f, .logInnov = 0.0f, .logW = 1.0f,
};
static wallRef_t rightWall = {
  .ref = 0.0f, .var = 0.01f, .lastMs = 0, .seeded = false,
  .posAxis = KC_STATE_Y, .rRow = 1, .sign = 1.0f, .dzSign = -1.0f,
  .logRef = 0.0f, .logInnov = 0.0f, .logW = 1.0f,
};

// Takeoff heading, captured from the quaternion while on the ground. The yaw-hold down-weight
// inflates R as the craft rotates away from this heading (the beam->axis mapping assumes it).
static float wallYaw0 = 0.0f;
static bool wallYaw0Seeded = false;

// diagnostics: the implied X/Y a single wall's reading places the craft at (ref - sign*r_h on the
// last-updated wall of each axis), for the shadow comparison against stateEstimate.x/y.
static float wallXImplied = 0.0f;
static float wallYImplied = 0.0f;

// Smallest representable angle difference yaw-yaw0 in [-pi, pi] (mirrors Python angle_diff).
static float wallAngleDiff(float a, float b)
{
  float d = a - b;
  while (d > (float)M_PI)  { d -= 2.0f * (float)M_PI; }
  while (d < -(float)M_PI) { d += 2.0f * (float)M_PI; }
  return d;
}

static void wallUpdate(kalmanCoreData_t* this, tofMeasurement_t* tof,
                       bool quadIsFlying, wallRef_t* wall)
{
  // Reject out-of-range / sentinel returns defensively (the deck driver also gates validity).
  if (!(tof->distance > 0.05f && tof->distance < 5.0f)) {
    return;
  }

  // --- per-beam geometry from the rotation matrix R (the world-Z of this body axis) ---
  //   dz is the world-Z component of the beam direction; cosSlant is the horizontal fraction
  //   sqrt(1 - dz^2). cosSlant clamped to >= wallRaycastCosMin (Python abs_raycast_cos_min), so a
  //   near-vertical beam (which the horizontal model cannot use) does not blow up r_h.
  float dz = wall->dzSign * this->R[2][wall->rRow];
  float cosSlant = 1.0f - dz * dz;
  if (cosSlant < 0.0f) { cosSlant = 0.0f; }
  cosSlant = sqrtf(cosSlant);
  if (cosSlant < wallRaycastCosMin) { cosSlant = wallRaycastCosMin; }
  float r_h = tof->distance * cosSlant;            // horizontal distance to the wall

  // --- reference random-walk predict (dt since THIS wall last updated; ~10 Hz per beam) -- a
  //     stiff datum that drifts only slowly, so a real wall does not wander. ---
  if (wall->lastMs != 0) {
    float dt = (float)(tof->timestamp - wall->lastMs) * 0.001f;
    if (dt > 0.0f) {
      if (dt > wallDtMaxS) { dt = wallDtMaxS; }
      wall->var += wallRwStd * wallRwStd * dt;
    }
  }
  wall->lastMs = tof->timestamp;

  // --- on the ground: re-capture the wall reference every frame (resilient to estimator reset /
  //     a prior flight's re-seat). The wall distance is unknown, so capture it from the beam
  //     (X,Y ~ 0 on the ground => wallRef ~ the preflight wall range). Also capture the takeoff
  //     heading yaw0 here (the beam->axis mapping is anchored to it). ---
  if (!quadIsFlying) {
    // measured_pos = wallRef + sign*r_h with pos == S[posAxis] => wallRef = S[posAxis] - sign*r_h
    wall->ref = this->S[wall->posAxis] - wall->sign * r_h;
    wall->var = wallInitVarM2;
    wall->seeded = true;
    wall->logRef = wall->ref;
    // yaw0 from the quaternion (atan2 form, mirrors kalman_core.c externalize)
    wallYaw0 = atan2f(2.0f * (this->q[1] * this->q[2] + this->q[0] * this->q[3]),
                      this->q[0] * this->q[0] + this->q[1] * this->q[1]
                        - this->q[2] * this->q[2] - this->q[3] * this->q[3]);
    wallYaw0Seeded = true;
    return;
  }
  if (!wall->seeded) {
    wall->ref = this->S[wall->posAxis] - wall->sign * r_h;
    wall->var = wallInitVarM2;
    wall->seeded = true;
  }

  // --- heteroscedastic measurement variance from the sensor's own per-shot sigma ---
  //   tof->stdDev carries the VL53L1x sigma [m]; floor it and add a model term in quadrature
  //   (matches Python measured_range_var). 0 / invalid -> a fixed fallback std.
  float rangeVar;
  if (tof->stdDev > 0.0f) {
    float sig = tof->stdDev;
    if (sig < wallSigFloorM) { sig = wallSigFloorM; }
    rangeVar = sig * sig + wallSigModelM * wallSigModelM;
  } else {
    rangeVar = wallStdFallbackM * wallStdFallbackM;
  }
  float measPosVar = rangeVar * cosSlant * cosSlant;  // project the range variance to the axis

  // --- YAW-HOLD / tilt down-weight: inflate R as the craft tilts or rotates from the takeoff
  //     heading. roll/pitch/yaw recovered from the quaternion (same atan2/asin form as
  //     kalman_core.c). yaw_dev is the angle from the captured yaw0. ---
  float yaw = atan2f(2.0f * (this->q[1] * this->q[2] + this->q[0] * this->q[3]),
                     this->q[0] * this->q[0] + this->q[1] * this->q[1]
                       - this->q[2] * this->q[2] - this->q[3] * this->q[3]);
  float pitch = asinf(-2.0f * (this->q[1] * this->q[3] - this->q[0] * this->q[2]));
  float roll = atan2f(2.0f * (this->q[2] * this->q[3] + this->q[0] * this->q[1]),
                      this->q[0] * this->q[0] - this->q[1] * this->q[1]
                        - this->q[2] * this->q[2] + this->q[3] * this->q[3]);
  float yawDev = wallYaw0Seeded ? wallAngleDiff(yaw, wallYaw0) : 0.0f;
  float worstDeg = fabsf(roll) * RAD_TO_DEG;
  float pitchDeg = fabsf(pitch) * RAD_TO_DEG;
  float yawDevDeg = fabsf(yawDev) * RAD_TO_DEG;
  if (pitchDeg > worstDeg)  { worstDeg = pitchDeg; }
  if (yawDevDeg > worstDeg) { worstDeg = yawDevDeg; }
  float scale = (wallTiltRScaleDeg > 1e-6f) ? wallTiltRScaleDeg : 1e-6f;
  float tiltMult = 1.0f + (worstDeg / scale) * (worstDeg / scale);

  // --- RAY-CAST wall-vs-floor down-weight: a beam that descends (dz < 0) reaches the floor
  //     (world Z=0) at slant s_floor = -z/dz; if that is nearer than the predicted wall slant the
  //     beam is probably on the floor, so inflate R by 1/w_wall (a sigmoid). The predicted wall
  //     slant MUST come from the STATE ((ref + sign*pos)/cosSlant, == Python sensors.py pred_wall) --
  //     NOT the measured range, which equals s_floor exactly when the beam IS on the floor and would
  //     then defeat the test (margin ~ 0 => no down-weight). dz >= 0 (level/up) => w_wall = 1. ---
  float floorRMult = 1.0f;
  if (dz < -1e-3f) {
    float z = this->S[KC_STATE_Z];
    float sFloor = -z / dz;                       // slant range at which the beam reaches Z=0
    float predWall = (wall->ref + wall->sign * this->S[wall->posAxis]) / cosSlant; // state-predicted wall slant
    float margin = (sFloor - predWall) / wallRaycastFloorScaleM;
    if (margin > 30.0f)  { margin = 30.0f; }
    if (margin < -30.0f) { margin = -30.0f; }
    float wWall = 1.0f / (1.0f + expf(-margin));  // ~1 wall clearly nearer, ~0 floor nearer
    float wFloor = (wWall > wallRaycastFloorWMin) ? wWall : wallRaycastFloorWMin;
    floorRMult = 1.0f / wFloor;
  }

  // --- position pseudo-measurement: measured_pos = ref + sign*r_h,  h[posAxis] = 1 ---
  float measuredPos = wall->ref + wall->sign * r_h;
  float innov = measuredPos - this->S[wall->posAxis];

  // Student-t robust weight. The innovation variance folds in BOTH the reference uncertainty and
  // the position uncertainty, so "how surprised are we" is measured against everything that is
  // genuinely uncertain (a freshly re-seated reference does not read as a furniture event). The
  // geometric defenses (tilt/yaw, ray-cast floor) inflate the measurement variance going in.
  float measPosVarEff = measPosVar * tiltMult * floorRMult;
  float Ppp = this->P[wall->posAxis][wall->posAxis];
  float Sinnov = measPosVarEff + wall->var + Ppp;
  float w = 1.0f;
  if (Sinnov > 0.0f) {
    float rstd2 = (innov * innov) / Sinnov;
    w = (wallNu + 1.0f) / (wallNu + rstd2);
  }
  float wClip = w;
  if (wClip > 1.0f) { wClip = 1.0f; }
  if (wClip < wallWMin) { wClip = wallWMin; }
  wall->logRef = wall->ref;
  wall->logInnov = innov;
  wall->logW = w;
  // the X/Y this wall alone places the craft at (for the shadow comparison) -- same form as
  // measuredPos above (ref + sign*r_h), NOT minus.
  float implied = wall->ref + wall->sign * r_h;
  if (wall->posAxis == KC_STATE_X) { wallXImplied = implied; }
  else                             { wallYImplied = implied; }

  // --- LAYER 1: main EKF position update, R inflated by 1/w (reference uncertainty included) ---
  float Reff = (measPosVarEff + wall->var) / wClip;
  float h[KC_STATE_DIM] = {0};
  arm_matrix_instance_f32 H = {1, KC_STATE_DIM, h};
  h[wall->posAxis] = 1.0f;
  kalmanCoreScalarUpdate(this, &H, innov, sqrtf(Reff));

  // --- LAYER 2: jump-diffusion re-seat of THIS wall's reference ---
  //   Inject variance only on a clear wall swap (w below reseatWMax) so a smooth translation
  //   (a moderate, persistent w) does NOT loosen the datum. Then a scalar KF step pulls the
  //   reference toward the position-implied observation. refObs uses the POST-update position: the
  //   reference re-seats only as far as the opposing wall / flow held the position.
  if (w < wallReseatWMax) {
    float inject = (1.0f - wClip) * wallJumpCouple * wallJumpSizeM * wallJumpSizeM;
    float room = wallAVarMaxM2 - wall->var;
    if (room < 0.0f) { room = 0.0f; }
    if (inject > room) { inject = room; }
    wall->var += inject;
  }
  // refObs observes the SAME quantity as the capture/seed (ref = S[posAxis] - sign*r_h), so it MUST
  // use minus -- a plus here re-seats the wall to the wrong side of the craft on every swap (runaway).
  float refObs = this->S[wall->posAxis] - wall->sign * r_h;
  float Kref = wall->var / (wall->var + measPosVarEff);
  wall->ref += Kref * (refObs - wall->ref);
  wall->var *= (1.0f - Kref);
}

// Clear all module-static wall state so the next flight re-seeds its references FRESH (at the first
// in-flight frame, with X,Y ~ 0 at takeoff) instead of inheriting a prior flight's re-seated box/hand
// reference. MUST be called on every estimator reset: unlike the down/up surfaces (whose down beam
// re-grounds every pre-flight frame), the horizontal beams do not reliably re-capture on the ground,
// so without this a second flight flies on stale references -- it drove into the table (flight
// 150720: fRef stuck at 0.35 from the previous hand-circle instead of the true 1.67).
void kalmanCoreWallReset(void)
{
  wallRef_t* walls[4] = { &frontWall, &backWall, &leftWall, &rightWall };
  for (int i = 0; i < 4; i++) {
    walls[i]->ref = 0.0f;
    walls[i]->var = 0.01f;   // init_ref_std_m^2 (matches the static initialisers above)
    walls[i]->lastMs = 0;
    walls[i]->seeded = false;
    walls[i]->logRef = 0.0f;
    walls[i]->logInnov = 0.0f;
    walls[i]->logW = 1.0f;
  }
  wallYaw0 = 0.0f;
  wallYaw0Seeded = false;
}

void kalmanCoreUpdateWithWallFront(kalmanCoreData_t* this, tofMeasurement_t* tof, bool quadIsFlying)
{
  wallUpdate(this, tof, quadIsFlying, &frontWall);
}

void kalmanCoreUpdateWithWallBack(kalmanCoreData_t* this, tofMeasurement_t* tof, bool quadIsFlying)
{
  wallUpdate(this, tof, quadIsFlying, &backWall);
}

void kalmanCoreUpdateWithWallLeft(kalmanCoreData_t* this, tofMeasurement_t* tof, bool quadIsFlying)
{
  wallUpdate(this, tof, quadIsFlying, &leftWall);
}

void kalmanCoreUpdateWithWallRight(kalmanCoreData_t* this, tofMeasurement_t* tof, bool quadIsFlying)
{
  wallUpdate(this, tof, quadIsFlying, &rightWall);
}

/**
 * Wall absolute-position (front/back/left/right beam) fusion tuning. Threshold-free: these shape
 * the Student-t robust weight, the self-calibrating reference and the geometric (tilt/yaw, floor
 * ray-cast) down-weights, NOT a gate. Defaults are ported from the validated Python wall config and
 * mirror the zsurf height defaults -- they should rarely need touching.
 */
PARAM_GROUP_START(zwall)
/**
 * @brief Student-t degrees of freedom: lower => heavier tails => SOFTER rejection of a large
 * innovation (a box pushed to a wall). The robust weight is w = (nu+1)/(nu + innov^2/S).
 */
PARAM_ADD(PARAM_FLOAT, nu, &wallNu)
/**
 * @brief Robust-weight floor: the most a beam's R can be inflated is 1/wMin. Smaller => a
 * disagreeing beam is more strongly (but never fully) ignored.
 */
PARAM_ADD(PARAM_FLOAT, wMin, &wallWMin)
/**
 * @brief Lower bound [m] on the per-shot VL53L1x sigma used for the measurement R.
 */
PARAM_ADD(PARAM_FLOAT, sigFloor, &wallSigFloorM)
/**
 * @brief Wall non-planarity / mount model std [m], added in quadrature to the sensor sigma.
 */
PARAM_ADD(PARAM_FLOAT, sigModel, &wallSigModelM)
/**
 * @brief Jump-diffusion coupling: fraction of jumpSize^2 injected into a wall's reference variance
 * as the robust weight -> 0 (a wall swap), letting the reference re-seat onto the new surface.
 */
PARAM_ADD(PARAM_FLOAT, jumpCouple, &wallJumpCouple)
/**
 * @brief Wall-swap step scale [m] for the jump-diffusion injection.
 */
PARAM_ADD(PARAM_FLOAT, jumpSize, &wallJumpSizeM)
/**
 * @brief Cap [m^2] on a reference's variance (bounds how fast it can re-seat per event).
 */
PARAM_ADD(PARAM_FLOAT, aVarMax, &wallAVarMaxM2)
/**
 * @brief Re-seat only when the robust weight is below this (a clear wall swap), so a smooth
 * translation (moderate persistent weight) does not loosen the datum.
 */
PARAM_ADD(PARAM_FLOAT, reseatW, &wallReseatWMax)
/**
 * @brief Per-wall reference random-walk std [m/sqrt(s)] -- how fast a wall datum may drift (a touch
 * faster than the stiff floor datum, walls move with the room layout).
 */
PARAM_ADD(PARAM_FLOAT, rw, &wallRwStd)
/**
 * @brief Yaw/tilt R-inflation scale [deg]: R *= 1 + (max(|roll|,|pitch|,|yawDev|)/this)^2, so the
 * walls fade out smoothly as the craft rotates away from the takeoff heading.
 */
PARAM_ADD(PARAM_FLOAT, tiltRScaleDeg, &wallTiltRScaleDeg)
/**
 * @brief Ray-cast wall/floor margin sigmoid scale [m]: how sharply a downward beam crosses over
 * from "hitting the wall" to "hitting the floor".
 */
PARAM_ADD(PARAM_FLOAT, floorScale, &wallRaycastFloorScaleM)
/**
 * @brief Clamp on the horizontal fraction cosSlant = sqrt(1 - dz^2): a near-vertical beam cannot
 * blow up the horizontal range r_h below this.
 */
PARAM_ADD(PARAM_FLOAT, cosMin, &wallRaycastCosMin)
/**
 * @brief Floor-confidence lower clamp: the most a floor-grazing beam's R can be inflated is
 * 1/floorWMin.
 */
PARAM_ADD(PARAM_FLOAT, floorWMin, &wallRaycastFloorWMin)
PARAM_GROUP_STOP(zwall)

/**
 * Wall absolute-position diagnostics: the self-calibrating per-wall references, the position
 * innovations, the robust weights, and the implied X/Y each axis's last wall placed the craft at
 * (the shadow comparison against stateEstimate.x/y during the OBSERVE flight). A weight collapsing
 * toward 0 with the reference then re-seating is the box-arriving re-association event; the opposing
 * wall / flow should hold the position.
 */
LOG_GROUP_START(zwall)
LOG_ADD(LOG_FLOAT, fRef, &frontWall.logRef)   // front-wall reference along X [m]
LOG_ADD(LOG_FLOAT, bRef, &backWall.logRef)    // back-wall reference along X [m]
LOG_ADD(LOG_FLOAT, lRef, &leftWall.logRef)    // left-wall reference along Y [m]
LOG_ADD(LOG_FLOAT, rRef, &rightWall.logRef)   // right-wall reference along Y [m]
LOG_ADD(LOG_FLOAT, fW, &frontWall.logW)       // front robust weight (tiny => swap)
LOG_ADD(LOG_FLOAT, bW, &backWall.logW)        // back robust weight
LOG_ADD(LOG_FLOAT, lW, &leftWall.logW)        // left robust weight
LOG_ADD(LOG_FLOAT, rW, &rightWall.logW)       // right robust weight
LOG_ADD(LOG_FLOAT, fInnov, &frontWall.logInnov) // front position innovation [m]
LOG_ADD(LOG_FLOAT, lInnov, &leftWall.logInnov)  // left position innovation [m]
LOG_ADD(LOG_FLOAT, xImplied, &wallXImplied)   // X the last-updated front/back wall implies [m]
LOG_ADD(LOG_FLOAT, yImplied, &wallYImplied)   // Y the last-updated left/right wall implies [m]
LOG_GROUP_STOP(zwall)
