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

#include "mm_tof_surface.h"
#include "param.h"
#include "log.h"
#include <math.h>

// ============================================================================
//  Opposing-surface absolute-height fusion (DOWN beam + UP beam)
//
//  Naming: the two references are named by the BEAM that sees them -- "down"
//  (zranger) and "up" (multiranger up) -- NOT "floor"/"ceiling". That is
//  deliberate: each reference SELF-CALIBRATES to whatever surface is actually in
//  the beam, so the down surface is usually the floor but becomes a TABLE when one
//  slides under, and the up surface is usually the ceiling but becomes a LANTERN /
//  a HAND. The beam direction is the only thing that is always true.
//
//  The one real asymmetry is the DATUM, not the surface kind: the down surface has
//  a known nominal (the ground at Z=0), while the up surface's height is unknown
//  and is captured on the ground at takeoff. That is what `captureOnGround` encodes.
//
//  A firmware port of the validated Python X/Y wall machinery
//  (fusion/sensors.py: robust_scalar_update + rangeabs_z_update) applied to the
//  vertical axis, adding NO Kalman states. Two ideas, both threshold-free -- they
//  replace the run-by-run hand-tuned gates that locked the truth out (e.g. flight
//  115242: a table re-referenced the height because its 0.50 m innovation exceeded
//  a 0.12 m gate, so the only correct sensor was rejected).
//
//   LAYER 1 -- Student-t robust scalar update (replaces every hard gate).
//     Each beam is a direct height pseudo-measurement measured_z = ref +/- vert.
//     Instead of "if |innov|>gate: reject", the innovation INFLATES R by 1/w via
//     the Student-t weight  w = (nu+1)/(nu + innov^2/S). A clean reading (small
//     innov) updates at full strength; furniture (large innov) is softly -- never
//     hard -- down-weighted, so the filter decides, not a cutoff. The per-shot
//     VL53L1x sigma feeds R directly (heteroscedastic): a grazing / weak / far
//     return self-reports high sigma and self-weakens.
//
//   LAYER 2 -- per-surface self-calibrating reference (auxiliary 1-D KF).
//     A surface height is a slow random-walk scalar (ref +/- var), not a frozen
//     capture. When a beam persistently disagrees with the height held by the
//     OPPOSING surface, the low robust weight injects variance into THAT surface's
//     reference (jump-diffusion) and the reference re-seats onto the new surface
//     (table -> down ref ~ table height; lantern -> up ref drops). Height never
//     follows the furniture -- the opposing surface anchors it through the
//     one-frame re-seat -- and the moved reference then makes the beam consistent
//     again, so the fast sensor is RECLAIMED, not merely ignored.
//
//  The re-seat observes the POST-update height, so a reference only moves to the
//  extent the height did NOT follow its beam (i.e. to the extent the opposing
//  surface anchored it). With a single beam (no opposing signal) height and
//  reference move together -- the irreducible single-beam ambiguity, handled
//  gracefully rather than papered over with a threshold.
//
//  Cascade caveat: running the reference KF sequentially with the main EKF (rather
//  than as joint states, which KC_STATE_DIM=9 has no room for) drops the shared
//  cross-covariance -- the standard Schmidt-style approximation. The Student-t R
//  inflation (1/wMin up to 1000x) and folding the reference variance into the
//  height-measurement R both suppress any height move during a swap, so z stays
//  put while the reference catches up.
// ============================================================================

// --- ported from the validated Python fusion config (fusion/state.py defaults) ---
static float surfNu          = 4.0f;    // Student-t dof (abs_dof): heavy tails -> soft rejection
static float surfWMin        = 0.001f;  // robust weight floor (abs_w_min): max R inflation = 1/wMin
static float surfSigFloorM   = 0.005f;  // meas_sigma_floor_mm: VL53L1x sigma is floored here [m]
static float surfSigModelM   = 0.020f;  // meas_sigma_model_mm: surface non-planarity / mount [m]
static float surfStdFallbackM = 0.020f; // abs_meas_std_m: used when the per-shot sigma is invalid
static float surfJumpCouple  = 0.5f;    // abs_jump_couple: fraction of jumpSize^2 injected at w->0
static float surfJumpSizeM   = 0.35f;   // abs_jump_size_m: surface-swap step scale [m]
static float surfAVarMaxM2   = 0.25f;   // abs_a_var_max_m2: cap on a reference's variance [m^2]
static float surfReseatWMax  = 0.40f;   // z_reseat_w_max: only re-seat below this w (a clear swap)
static float surfDownRwStd   = 0.005f;  // down-surface reference random-walk std [m/sqrt(s)]
static float surfUpRwStd     = 0.005f;  // up-surface reference random-walk std [m/sqrt(s)]
static float surfInitVarM2   = 0.01f;   // init_ref_std_m^2: reference variance at (re)seed [m^2]
static const float surfDtMaxS = 0.20f;  // clamp the random-walk dt across a sensor gap [s]

typedef struct {
  float ref;            // surface height above ground [m] (down ~0; up ~room height)
  float var;            // reference variance [m^2]
  uint32_t lastMs;      // tick of the last update (random-walk dt)
  bool seeded;          // reference initialised
  bool captureOnGround; // datum asymmetry: up captures its height from the beam on the ground;
                        // down's nominal IS the ground (Z=0, known) so it is not captured.
  float sign;           // +1 down (z = ref + vert), -1 up (z = ref - vert)
  // diagnostics (first flights)
  float logRef;         // current reference [m]
  float logInnov;       // height innovation [m]
  float logW;           // robust weight (<1 => down-weighted; tiny => surface swap)
} surfaceRef_t;

// Module-static reference state (formerly kalmanCoreData_t::ceilingRef; moved here so a surface
// reference is local to its fusion and re-grounds itself on the ground -- see below). Not reset by
// kalmanCoreInit, but the on-ground re-grounding below makes that unnecessary: every pre-flight
// frame forces the down ref to 0 and re-captures the up ref, so a prior flight's re-seated
// table/lantern reference cannot persist across an estimator reset.
static surfaceRef_t downSurf = {
  .ref = 0.0f, .var = 0.01f, .lastMs = 0, .seeded = true,
  .captureOnGround = false, .sign = 1.0f, .logRef = 0.0f, .logInnov = 0.0f, .logW = 1.0f,
};
static surfaceRef_t upSurf = {
  .ref = 0.0f, .var = 0.01f, .lastMs = 0, .seeded = false,
  .captureOnGround = true, .sign = -1.0f, .logRef = 0.0f, .logInnov = 0.0f, .logW = 1.0f,
};

static void surfaceUpdate(kalmanCoreData_t* this, tofMeasurement_t* tof,
                          bool quadIsFlying, surfaceRef_t* surf)
{
  // Skip near-horizontal (the height model diverges as R[2][2] -> 0), exactly as mm_tof.
  if (!(fabsf(this->R[2][2]) > 0.1f && this->R[2][2] > 0.0f)) {
    return;
  }
  // Reject out-of-range / sentinel returns defensively (the deck driver also gates validity).
  if (!(tof->distance > 0.05f && tof->distance < 5.0f)) {
    return;
  }

  float angle = fabsf(acosf(this->R[2][2])) - DEG_TO_RAD * (15.0f / 2.0f);
  if (angle < 0.0f) {
    angle = 0.0f;
  }
  float cosAngle = cosf(angle);
  float vert = tof->distance * cosAngle;          // vertical component of the slant range
  float rwStd = surf->captureOnGround ? surfUpRwStd : surfDownRwStd;

  // --- reference random-walk predict (dt since THIS surface last updated; down ~40 Hz, up
  //     ~10 Hz) -- a stiff datum that drifts only slowly, so a real surface does not wander. ---
  if (surf->lastMs != 0) {
    float dt = (float)(tof->timestamp - surf->lastMs) * 0.001f;
    if (dt > 0.0f) {
      if (dt > surfDtMaxS) { dt = surfDtMaxS; }
      surf->var += rwStd * rwStd * dt;
    }
  }
  surf->lastMs = tof->timestamp;

  // --- on the ground: re-ground the reference every frame (resilient to estimator reset / a prior
  //     flight's re-seat). DOWN: the surface IS the ground, ref == 0 by definition. UP: the
  //     room/lantern height is unknown, so capture it from the beam (Z ~ 0 on the ground). ---
  if (!quadIsFlying) {
    surf->ref = surf->captureOnGround ? (this->S[KC_STATE_Z] - surf->sign * vert) : 0.0f;
    surf->var = surfInitVarM2;
    surf->seeded = true;
    surf->logRef = surf->ref;
    return;
  }
  if (!surf->seeded) {
    surf->ref = this->S[KC_STATE_Z] - surf->sign * vert;
    surf->var = surfInitVarM2;
    surf->seeded = true;
  }

  // --- heteroscedastic measurement variance from the sensor's own per-shot sigma ---
  //   tof->stdDev carries the VL53L1x sigma [m]; floor it and add a model term in quadrature
  //   (matches Python measured_range_var). 0 / invalid -> a fixed fallback std.
  float rangeVar;
  if (tof->stdDev > 0.0f) {
    float sig = tof->stdDev;
    if (sig < surfSigFloorM) { sig = surfSigFloorM; }
    rangeVar = sig * sig + surfSigModelM * surfSigModelM;
  } else {
    rangeVar = surfStdFallbackM * surfStdFallbackM;
  }
  float measZVar = rangeVar * cosAngle * cosAngle;  // project the range variance to the vertical

  // --- height pseudo-measurement: measured_z = ref +/- vert,  h[Z] = 1 ---
  float measuredZ = surf->ref + surf->sign * vert;
  float innov = measuredZ - this->S[KC_STATE_Z];

  // Student-t robust weight. The innovation variance folds in BOTH the reference uncertainty and
  // the height uncertainty, so "how surprised are we" is measured against everything that is
  // genuinely uncertain (a freshly re-seated reference does not read as a furniture event).
  float Pzz = this->P[KC_STATE_Z][KC_STATE_Z];
  float Sinnov = measZVar + surf->var + Pzz;
  float w = 1.0f;
  if (Sinnov > 0.0f) {
    float rstd2 = (innov * innov) / Sinnov;
    w = (surfNu + 1.0f) / (surfNu + rstd2);
  }
  float wClip = w;
  if (wClip > 1.0f) { wClip = 1.0f; }
  if (wClip < surfWMin) { wClip = surfWMin; }
  surf->logRef = surf->ref;
  surf->logInnov = innov;
  surf->logW = w;

  // --- LAYER 1: main EKF height update, R inflated by 1/w (reference uncertainty included) ---
  float Reff = (measZVar + surf->var) / wClip;
  float h[KC_STATE_DIM] = {0};
  arm_matrix_instance_f32 H = {1, KC_STATE_DIM, h};
  h[KC_STATE_Z] = 1.0f;
  kalmanCoreScalarUpdate(this, &H, innov, sqrtf(Reff));

  // --- LAYER 2: jump-diffusion re-seat of THIS surface's reference ---
  //   Inject variance only on a clear surface swap (w below reseatWMax) so a smooth climb/descent
  //   (a moderate, persistent w ~ 0.7) does NOT loosen the datum. Then a scalar KF step pulls the
  //   reference toward the height-implied observation. refObs uses the POST-update height: the
  //   reference re-seats only as far as the opposing surface held the height.
  if (w < surfReseatWMax) {
    float inject = (1.0f - wClip) * surfJumpCouple * surfJumpSizeM * surfJumpSizeM;
    float room = surfAVarMaxM2 - surf->var;
    if (room < 0.0f) { room = 0.0f; }
    if (inject > room) { inject = room; }
    surf->var += inject;
  }
  float refObs = this->S[KC_STATE_Z] - surf->sign * vert;
  float Kref = surf->var / (surf->var + measZVar);
  surf->ref += Kref * (refObs - surf->ref);
  surf->var *= (1.0f - Kref);
}

void kalmanCoreUpdateWithTofDown(kalmanCoreData_t* this, tofMeasurement_t* tof, bool quadIsFlying)
{
  surfaceUpdate(this, tof, quadIsFlying, &downSurf);
}

void kalmanCoreUpdateWithTofUp(kalmanCoreData_t* this, tofMeasurement_t* tof, bool quadIsFlying)
{
  surfaceUpdate(this, tof, quadIsFlying, &upSurf);
}

/**
 * Opposing-surface (down + up beam) height-fusion tuning. Threshold-free: these shape the
 * Student-t robust weight and the self-calibrating reference, NOT a gate. Defaults are ported
 * from the validated Python wall config -- they should rarely need touching.
 */
PARAM_GROUP_START(zsurf)
/**
 * @brief Student-t degrees of freedom: lower => heavier tails => SOFTER rejection of a large
 * innovation (furniture). The robust weight is w = (nu+1)/(nu + innov^2/S).
 */
PARAM_ADD(PARAM_FLOAT, nu, &surfNu)
/**
 * @brief Robust-weight floor: the most a beam's R can be inflated is 1/wMin. Smaller => a
 * disagreeing beam is more strongly (but never fully) ignored.
 */
PARAM_ADD(PARAM_FLOAT, wMin, &surfWMin)
/**
 * @brief Lower bound [m] on the per-shot VL53L1x sigma used for the measurement R.
 */
PARAM_ADD(PARAM_FLOAT, sigFloor, &surfSigFloorM)
/**
 * @brief Surface non-planarity / mount model std [m], added in quadrature to the sensor sigma.
 */
PARAM_ADD(PARAM_FLOAT, sigModel, &surfSigModelM)
/**
 * @brief Jump-diffusion coupling: fraction of jumpSize^2 injected into a surface's reference
 * variance as the robust weight -> 0 (a surface swap), letting the reference re-seat.
 */
PARAM_ADD(PARAM_FLOAT, jumpCouple, &surfJumpCouple)
/**
 * @brief Surface-swap step scale [m] for the jump-diffusion injection.
 */
PARAM_ADD(PARAM_FLOAT, jumpSize, &surfJumpSizeM)
/**
 * @brief Cap [m^2] on a reference's variance (bounds how fast it can re-seat per event).
 */
PARAM_ADD(PARAM_FLOAT, aVarMax, &surfAVarMaxM2)
/**
 * @brief Re-seat only when the robust weight is below this (a clear surface swap), so a smooth
 * climb/descent (moderate persistent weight) does not loosen the datum.
 */
PARAM_ADD(PARAM_FLOAT, reseatW, &surfReseatWMax)
/**
 * @brief Down-surface reference random-walk std [m/sqrt(s)] -- how fast the down datum may drift.
 */
PARAM_ADD(PARAM_FLOAT, downRw, &surfDownRwStd)
/**
 * @brief Up-surface reference random-walk std [m/sqrt(s)] -- how fast the up datum may drift.
 */
PARAM_ADD(PARAM_FLOAT, upRw, &surfUpRwStd)
PARAM_GROUP_STOP(zsurf)

/**
 * Opposing-surface height diagnostics: the self-calibrating references, the height innovations,
 * and the robust weights. A weight collapsing toward 0 with the reference then re-seating is the
 * furniture (table/lantern/hand) re-association event; the opposing surface should hold the height.
 */
LOG_GROUP_START(zsurf)
LOG_ADD(LOG_FLOAT, dRef, &downSurf.logRef)     // down-surface reference [m] (~0; jumps to a table)
LOG_ADD(LOG_FLOAT, uRef, &upSurf.logRef)       // up-surface reference [m] (room; drops to a lantern)
LOG_ADD(LOG_FLOAT, dInnov, &downSurf.logInnov) // down height innovation [m]
LOG_ADD(LOG_FLOAT, uInnov, &upSurf.logInnov)   // up height innovation [m]
LOG_ADD(LOG_FLOAT, dW, &downSurf.logW)         // down robust weight (tiny => swap)
LOG_ADD(LOG_FLOAT, uW, &upSurf.logW)           // up robust weight (tiny => swap)
LOG_GROUP_STOP(zsurf)
