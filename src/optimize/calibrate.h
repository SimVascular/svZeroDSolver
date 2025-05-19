// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file calibrate.h
 * @brief opt::calibrate source file
 */

#ifndef SVZERODSOLVER_OPTIMIZE_CALIBRATOR_HPP_
#define SVZERODSOLVER_OPTIMIZE_CALIBRATOR_HPP_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <fstream>
#include <nlohmann/json.hpp>

#include "Model.h"
#include "debug.h"

/**
 * @brief Main function to run the 0D model calibration.
 * @param config JSON configuration for 0D model
 * @return Calibrated JSON configuration for the 0D model
 */
nlohmann::json calibrate(const nlohmann::json &config);

#endif  // SVZERODSOLVER_OPTIMIZE_CALIBRATOR_HPP_
