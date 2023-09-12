
#include "Parameter.h"

namespace zd_model {

Parameter::Parameter(int id, double value) 
{
  this->id = id;
  update(value);
}

Parameter::Parameter(int id, const std::vector<double> &times,
    const std::vector<double> &values, bool periodic) 
{
  this->id = id;
  this->isperiodic = periodic;
  update(times, values);
}

void Parameter::update(double value) 
{
  this->isconstant = true;
  this->isperiodic = true;
  this->value = value;
}


void Parameter::update(const std::vector<double> &times, const std::vector<double> &values) {
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

Parameter::~Parameter() {}

double Parameter::get(double time) 
{
  // Return the constant value if parameter is constant
  if (isconstant) {
    return value;
  }

  // Determine the time within this->times (necessary to extrapolate)
  double rtime;

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

void Parameter::to_steady() 
{
  if (isconstant) {
    return;
  }

  value = std::accumulate(values.begin(), values.end(), 0.0) / double(size);
  isconstant = true;
  steady_converted = true;
}

void Parameter::to_unsteady() 
{
  if (steady_converted) {
    isconstant = false;
    steady_converted = false;
  }
}

} 

