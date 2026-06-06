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

#pragma once

#include "kalman_core.h"

// Wall absolute-position fusion: the four horizontal multiranger beams (front/back/left/right)
// are fused as direct world-position pseudo-measurements with a Student-t robust weight (no hard
// gate) and a per-wall self-calibrating reference (an auxiliary 1-D KF that re-seats onto a box /
// pushed furniture via jump-diffusion). A faithful 4-wall mirror of mm_tof_surface.c: front/back
// anchor world X, left/right anchor world Y, with the per-beam sign and the wall reference standing
// in for the surface height. Adds NO Kalman states (KC_STATE_DIM stays 9) -- the references are
// module-static, exactly like the down/up surfaces.
//
// A firmware port of the validated Python wall machinery (fusion/sensors.py: rangeabs_update).
// tof->distance carries the slant range [m]; tof->stdDev carries the per-shot VL53L1x sigma [m]
// (0 => fall back to a fixed std); quadIsFlying drives the on-ground wall capture (each wall ref is
// captured from its beam at takeoff, X,Y ~ 0) and the yaw0 takeoff-heading capture.
// Clear all module-static wall references + seeded flags. MUST be called on every estimator reset so
// each flight re-seeds fresh -- the references persist across flights and are otherwise stale (a
// prior flight's re-seated box/hand reference would fly the next flight into a wall).
void kalmanCoreWallReset(void);

void kalmanCoreUpdateWithWallFront(kalmanCoreData_t* this, tofMeasurement_t* tof, bool quadIsFlying);
void kalmanCoreUpdateWithWallBack(kalmanCoreData_t* this, tofMeasurement_t* tof, bool quadIsFlying);
void kalmanCoreUpdateWithWallLeft(kalmanCoreData_t* this, tofMeasurement_t* tof, bool quadIsFlying);
void kalmanCoreUpdateWithWallRight(kalmanCoreData_t* this, tofMeasurement_t* tof, bool quadIsFlying);
