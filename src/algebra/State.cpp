// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "State.h"

State::State() {}

State::State(int n) {
  y = Eigen::Matrix<double, Eigen::Dynamic, 1>(n);
  ydot = Eigen::Matrix<double, Eigen::Dynamic, 1>(n);
}

State::~State() {}

State::State(const State &state) {
  y = state.y;
  ydot = state.ydot;
}

State State::Zero(int n) {
  // [TODO] what's going on here, returing a static State?
  static State state(n);
  state.y = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(n);
  state.ydot = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(n);
  return state;
}
