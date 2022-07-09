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
  TimeDependentParameter(std::vector<T> times, std::vector<T> values);

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
   * @brief Get the parameter value at the specified time.
   *
   * @param time Current time
   * @return Value at the time
   */
  T get(T time);

  bool isconstant;

  void to_steady();  ///< Bool value indicating if the parameter is constant
};

template <typename T>
TimeDependentParameter<T>::TimeDependentParameter() {}

template <typename T>
TimeDependentParameter<T>::TimeDependentParameter(std::vector<T> times,
                                                  std::vector<T> values) {
  this->times = times;
  this->values = values;
  size = times.size();
  if (size == 1) {
    isconstant = true;
    cycle_period = 1.0;
  } else {
    isconstant = false;
    cycle_period = times.back() - times[0];
  }
}

template <typename T>
TimeDependentParameter<T>::~TimeDependentParameter() {}

template <typename T>
T TimeDependentParameter<T>::get(T time) {
  // Return the first and only value if parameter is constant
  if (size == 1) {
    return values[0];
  }

  // Determine the time within a cycle (necessary to extrapolate)
  T rtime = fmod(time, cycle_period);

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
  values = std::vector<T>();
  times = std::vector<T>();
  values.push_back(mean);
  size = 1;
  isconstant = true;
}
}  // namespace MODEL

#endif  // SVZERODSOLVER_MODEL_PARAMETER_HPP_