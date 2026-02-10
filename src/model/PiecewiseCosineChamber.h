// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

/**
 * @file PiecewiseCosineChamber.h
 * @brief Backward compatibility alias for LinearElastanceChamber
 */

#ifndef SVZERODSOLVER_MODEL_PIECEWISECOSINECHAMBER_HPP_
#define SVZERODSOLVER_MODEL_PIECEWISECOSINECHAMBER_HPP_

#include "LinearElastanceChamber.h"

/**
 * @brief Backward compatibility alias for LinearElastanceChamber
 * 
 * PiecewiseCosineChamber has been renamed to LinearElastanceChamber.
 * This alias is maintained for backward compatibility with existing configurations.
 */
using PiecewiseCosineChamber = LinearElastanceChamber;

#endif  // SVZERODSOLVER_MODEL_PIECEWISECOSINECHAMBER_HPP_
