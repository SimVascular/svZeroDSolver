// Copyright (c) Stanford University, The Regents of the University of
//               California, and others.
//
// All Rights Reserved.
//
// See Copyright-SimVascular.txt for additional details.
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject
// to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
// OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/**
 * @file Parameter.h
 * @brief model::Parameter source file
 */
#ifndef SVZERODSOLVER_MODEL_PARAMETER_HPP_
#define SVZERODSOLVER_MODEL_PARAMETER_HPP_

#include <math.h>

#include <numeric>
#include <vector>

#include "DOFHandler.h"

/**
 * @brief Model Parameter.
 *
 * This class handles constant parameters and time-dependent parameters that
 * need to be interpolated and periodically applied.
 *
 */
class Parameter {
 public:
  /**
   * @brief Construct a new Parameter object
   *
   * @param id Global ID of the parameter
   * @param value The value of the parameter
   */
  Parameter(int id, double value);

  /**
   * @brief Construct a new Parameter object
   *
   * @param id Global ID of the parameter
   * @param times Time steps corresponding to the time-dependent values
   * @param values Values corresponding to the time steps
   */
  Parameter(int id, const std::vector<double>& times,
            const std::vector<double>& values, bool periodic = true);

  /**
   * @brief Destroy the Parameter object
   *
   */
  ~Parameter();

  int id;                      ///< Global ID of the parameter
  std::vector<double> times;   ///< Time steps if parameter is time-dependent
  std::vector<double> values;  ///< Values if parameter is time-dependent
  double value;                ///< Value if parameter is constant
  double cycle_period;  ///< Cardiac cycle period corresponding to the time
                        ///< sequence
  int size;         ///< Size of the time series if parameter is time-dependent
  bool isconstant;  ///< Bool value indicating if the parameter is constant
  bool isperiodic;  ///< Bool value indicating if the parameter is periodic with
                    ///< the cardiac cycle

  /**
   * @brief Update the parameter
   *
   * @param value Value of the parameter
   */
  void update(double value);

  /**
   * @brief Update the parameter
   *
   * @param times Time steps corresponding to the values
   * @param values Values correspondong to the time steps
   */
  void update(const std::vector<double>& times,
              const std::vector<double>& values);

  /**
   * @brief Get the parameter value at the specified time.
   *
   * @param time Current time
   * @return Value at the time
   */
  double get(double time);

  /**
   * @brief Convert the parameter into a steady mean state.
   *
   */
  void to_steady();

  /**
   * @brief Convert the parameter back into an unsteady mean state.
   *
   */
  void to_unsteady();

 private:
  bool steady_converted = false;
};

#endif  // SVZERODSOLVER_MODEL_PARAMETER_HPP_
