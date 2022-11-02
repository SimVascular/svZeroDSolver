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
 * @file timedependentparameter.hpp
 * @brief MODEL::TimeDependentParameter source file
 */
#ifndef SVZERODSOLVER_MODEL_PARAMETER_HPP_
#define SVZERODSOLVER_MODEL_PARAMETER_HPP_

#include <math.h>

#include <numeric>
#include <vector>

#include "dofhandler.hpp"

namespace MODEL {

/**
 * @brief Time-dependent parameter.
 *
 * This class handles time-dependent parameter that need to be interpolated
 * and periodically applied.
 *
 * @tparam T Scalar type (e.g. `float`, `double`)
 */
template <typename T>
class TimeDependentParameter {
 public:
  /**
   * @brief Construct a new Time Dependent Parameter object
   *
   */
  TimeDependentParameter();

  /**
   * @brief Construct a new Time Dependent Parameter object
   *
   * @param times Time steps corresponding to the values
   * @param values Values correspondong to the time steps
   */
  TimeDependentParameter(std::vector<T> times, std::vector<T> values, bool periodic = true);

  /**
   * @brief Destroy the Time Dependent Parameter object
   *
   */
  ~TimeDependentParameter();

  std::vector<T> times;   ///< Time steps
  std::vector<T> values;  ///< Values
  T cycle_period;  ///< Cardiac cycle period corresponding to the time sequence
  int size;

  /**
   * @brief Update parameters
   *
   * @param times Time steps corresponding to the values
   * @param values Values correspondong to the time steps
   */

  void update_params(std::vector<T> times, std::vector<T> values);

  /**
   * @brief Get the parameter value at the specified time.
   *
   * @param time Current time
   * @return Value at the time
   */
  T get(T time);

  bool isconstant;  ///< Bool value indicating if the parameter is constant
  bool isperiodic;  ///< Bool value indicating if the parameter is periodic with the cardiac cycle

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
  std::vector<T> times_cache;   ///< Time steps cache
  std::vector<T> values_cache;  ///< Values cache
  int size_cache;               ///< Size cache
};

template <typename T>
TimeDependentParameter<T>::TimeDependentParameter() {}

template <typename T>
TimeDependentParameter<T>::TimeDependentParameter(std::vector<T> times, std::vector<T> values, bool periodic) {
  this->times = times;
  this->values = values;
  size = values.size();
  if (size == 1) {
    isconstant = true;
    cycle_period = 1.0;
    isperiodic = true;
  } else {
    isconstant = false;
    cycle_period = times.back() - times[0];
    isperiodic = periodic;
  }
}

template <typename T>
void TimeDependentParameter<T>::update_params(std::vector<T> times, std::vector<T> values) {
  this->times = times;
  this->values = values;
  size = values.size();
  cycle_period = times.back() - times[0];
}

template <typename T>
TimeDependentParameter<T>::~TimeDependentParameter() {}

template <typename T>
T TimeDependentParameter<T>::get(T time) {
  // Return the first and only value if parameter is constant
  if (isconstant) {
    return values[0];
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
  return values[l] +
         ((values[k] - values[l]) / (times[k] - times[l])) * (rtime - times[l]);
}

template <typename T>
void TimeDependentParameter<T>::to_steady() {
  T mean =
      std::accumulate(values.begin(), values.end(), 0.0) / T(values.size());
  values_cache = values;
  times_cache = times;
  size_cache = size;
  values = std::vector<T>();
  times = std::vector<T>();
  values.push_back(mean);
  values.push_back(0.0);
  size = 1;
  isconstant = true;
}

template <typename T>
void TimeDependentParameter<T>::to_unsteady() {
  values = values_cache;
  times = times_cache;
  size = size_cache;
  if (size > 1) {
    isconstant = false;
  }
}

}  // namespace MODEL

#endif  // SVZERODSOLVER_MODEL_PARAMETER_HPP_
