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

#include "Parameter.h"

Parameter::Parameter(int id, double value) {
  this->id = id;
  update(value);
}

Parameter::Parameter(int id, const std::vector<double> &times,
                     const std::vector<double> &values, bool periodic) {
  this->id = id;
  this->isperiodic = periodic;
  update(times, values);
}

void Parameter::update(double value) {
  this->isconstant = true;
  this->isperiodic = true;
  this->value = value;
}

void Parameter::update(const std::vector<double> &times,
                       const std::vector<double> &values) {
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

double Parameter::get(double time) {
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

void Parameter::to_steady() {
  if (isconstant) {
    return;
  }

  value = std::accumulate(values.begin(), values.end(), 0.0) / double(size);
  isconstant = true;
  steady_converted = true;
}

void Parameter::to_unsteady() {
  if (steady_converted) {
    isconstant = false;
    steady_converted = false;
  }
}

