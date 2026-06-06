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

// Opposing-surface absolute-height fusion: the DOWN beam (zranger) and the UP beam (multiranger
// up) are fused as direct height pseudo-measurements with a Student-t robust weight (no hard gate)
// and a per-surface self-calibrating reference (an auxiliary 1-D KF that re-seats onto furniture
// via jump-diffusion). References are named by BEAM, not "floor"/"ceiling": each tracks whatever
// surface is in the beam (a table re-seats the down ref; a lantern / hand re-seats the up ref).
// A firmware port of the validated Python wall machinery (fusion/sensors.py), adding NO Kalman
// states. tof->distance carries the slant range [m]; tof->stdDev carries the per-shot VL53L1x
// sigma [m] (0 => fall back to a fixed std); quadIsFlying drives the on-ground reference
// (re-)grounding (down ref -> 0) / up-reference capture.
void kalmanCoreUpdateWithTofDown(kalmanCoreData_t* this, tofMeasurement_t* tof, bool quadIsFlying);
void kalmanCoreUpdateWithTofUp(kalmanCoreData_t* this, tofMeasurement_t* tof, bool quadIsFlying);
