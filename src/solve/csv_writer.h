// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file csv_writer.h
 * @brief csv_writer source file
 */
#ifndef SVZERODSOLVER_IO_CSVWRITER_HPP_
#define SVZERODSOLVER_IO_CSVWRITER_HPP_

#include <fstream>
#include <string>
#include <vector>

#include "Model.h"
#include "State.h"

std::string to_variable_csv(const std::vector<double> &times,
                            const std::vector<State> &states,
                            const Model &model, bool mean = false,
                            bool derivative = false);

std::string to_vessel_csv(const std::vector<double> &times,
                          const std::vector<State> &states, const Model &model,
                          bool mean = false, bool derivative = false);

#endif  // SVZERODSOLVER_IO_CSVWRITER_HPP_
