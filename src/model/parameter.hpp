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
 * @file parameter.hpp
 * @brief MODEL::Parameter source file
 */
#ifndef SVZERODSOLVER_MODEL_PARAMETER_HPP_
#define SVZERODSOLVER_MODEL_PARAMETER_HPP_

#include <math.h>

#include <numeric>
#include <vector>

#include "dofhandler.hpp"

namespace MODEL {

/**
 * @brief Model Parameter.
 *
 * This class handles constant parameters and time-dependent parameters that
 * need to be interpolated and periodically applied.
 *
 * @tparam T Scalar type (e.g. `float`, `double`)
 */
template <typename T>
class Parameter {
 public:
  /**
   * @brief Construct a new Parameter object
   *
   */
  Parameter();

  /**
   * @brief Construct a new Parameter object
   *
   * @param value The value of the parameter
   */
  Parameter(T value);

  /**
   * @brief Construct a new Parameter object
   *
   * @param times Time steps corresponding to the time-dependent values
   * @param values Values corresponding to the time steps
   */
  Parameter(std::vector<T> times, std::vector<T> values, bool periodic = true);

  /**
   * @brief Destroy the Parameter object
   *
   */
  ~Parameter();

  std::vector<T> times;   ///< Time steps if parameter is time-dependent
  std::vector<T> values;  ///< Values if parameter is time-dependent
  T value;                ///< Value if parameter is constant
  T cycle_period;   ///< Cardiac cycle period corresponding to the time sequence
  int size;         ///< Size of the time series if parameter is time-dependent
  bool isconstant;  ///< Bool value indicating if the parameter is constant
  bool isperiodic;  ///< Bool value indicating if the parameter is periodic with
                    ///< the cardiac cycle

  /**
   * @brief Update the parameter
   *
   * @param value Value of the parameter
   */
  void update(T value);

  /**
   * @brief Update the parameter
   *
   * @param times Time steps corresponding to the values
   * @param values Values correspondong to the time steps
   */

  void update(std::vector<T> times, std::vector<T> values);

  /**
   * @brief Get the parameter value at the specified time.
   *
   * @param time Current time
   * @return Value at the time
   */
  T get(T time);

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

template <typename T>
Parameter<T>::Parameter() {}

template <typename T>
Parameter<T>::Parameter(T value) {
  update(value);
}

template <typename T>
Parameter<T>::Parameter(std::vector<T> times, std::vector<T> values,
                        bool periodic) {
  this->isperiodic = periodic;
  update(times, values);
}

template <typename T>
void Parameter<T>::update(T value) {
  this->isconstant = true;
  this->isperiodic = true;
  this->value = value;
}

template <typename T>
void Parameter<T>::update(std::vector<T> times, std::vector<T> values) {
  this->size = values.size();
  if (this->size == 1) {
    this->value = values[0];
    this->isconstant = true;
  } else {
    this->times = times;
    this->values = values;
    this->cycle_period = times.back() - times[0];
    this->isconstant = false;
  }
}

template <typename T>
Parameter<T>::~Parameter() {}

template <typename T>
T Parameter<T>::get(T time) {
  // Return the constant value if parameter is constant
  if (isconstant) {
    return value;
  }

  // Determine the time within this->times (necessary to extrapolate)
  T rtime;
  if (isperiodic == true) {
    rtime = fmod(time, cycle_period);
  } else {
    // this->times is not periodic when running with external solver
    rtime = time;
  }

  // Determine the lower and upper element for interpolation
  auto i = lower_bound(times.begin(), times.end(), rtime);
  unsigned int k = i - times.begin();
  if (i == times.end())
    --i;
  else if (*i == rtime) {
    return values[k];
  }
  unsigned int l = k ? k - 1 : 1;

  // Perform linear interpolation
  // TODO: Implement periodic cubic spline
  return values[l] +
         ((values[k] - values[l]) / (times[k] - times[l])) * (rtime - times[l]);
}

template <typename T>
void Parameter<T>::to_steady() {
  if (isconstant) {
    return;
  }
  value = std::accumulate(values.begin(), values.end(), 0.0) / T(size);
  ;
  isconstant = true;
  steady_converted = true;
}

template <typename T>
void Parameter<T>::to_unsteady() {
  if (steady_converted) {
    isconstant = false;
    steady_converted = false;
  }
}

}  // namespace MODEL

#endif  // SVZERODSOLVER_MODEL_PARAMETER_HPP_
