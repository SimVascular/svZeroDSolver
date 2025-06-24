// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file BlockType.h
 * @brief Specifies the types of blocks and their parameters
 */
#ifndef SVZERODSOLVER_MODEL_BLOCK_TYPE_HPP_
#define SVZERODSOLVER_MODEL_BLOCK_TYPE_HPP_

#include <string>

/**
 * @brief The types of blocks supported by the solver
 */
enum class BlockType {
  blood_vessel = 0,
  junction = 1,
  blood_vessel_junction = 2,
  resistive_junction = 3,
  flow_bc = 4,
  pressure_bc = 5,
  resistance_bc = 6,
  windkessel_bc = 7,
  open_loop_coronary_bc = 8,
  closed_loop_coronary_left_bc = 9,
  closed_loop_coronary_right_bc = 10,
  closed_loop_rcr_bc = 11,
  closed_loop_heart_pulmonary = 12,
  valve_tanh = 13,
  chamber_elastance_inductor = 14
};

/**
 * @brief The classes/categories of blocks supported. Some classes require
 * special handling (e.g. closed_loop).
 */
enum class BlockClass {
  vessel = 0,
  junction = 1,
  boundary_condition = 2,
  closed_loop = 3,
  external = 4,
  valve = 5,
  chamber = 6
};

/**
 * @brief The types of vessel blocks supported.
 */
enum class VesselType { inlet = 0, outlet = 1, both = 2, neither = 3 };

#endif
