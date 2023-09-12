#include "State.h"

namespace algebra {

State::State() {}

State::State(unsigned int n) 
{
  y = Eigen::Matrix<double, Eigen::Dynamic, 1>(n);
  ydot = Eigen::Matrix<double, Eigen::Dynamic, 1>(n);
}

State::~State() {}

State::State(const State &state) 
{
  y = state.y;
  ydot = state.ydot;
}

State State::Zero(unsigned int n) 
{
  static State state(n);
  state.y = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(n);
  state.ydot = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(n);
  return state;
}

}  

