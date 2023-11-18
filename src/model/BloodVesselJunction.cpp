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

#include "BloodVesselJunction.h"

void BloodVesselJunction::setup_dofs(DOFHandler &dofhandler) {
  if (inlet_nodes.size() != 1) {
    throw std::runtime_error(
        "Blood vessel junction does not support multiple inlets.");
  }

  num_outlets = outlet_nodes.size();
  Block::setup_dofs_(dofhandler, num_outlets + 1, {});
  num_triplets.F = 1 + 4 * num_outlets;
  num_triplets.E = 3 * num_outlets;
  num_triplets.D = 2 * num_outlets;
}

void BloodVesselJunction::update_constant(SparseSystem &system,
                                          std::vector<double> &parameters) {
  // Mass conservation
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[1]) = 1.0;

  for (size_t i = 0; i < num_outlets; i++) {
    double inductance = parameters[global_param_ids[num_outlets + i]];
    double resistance = parameters[global_param_ids[i]];
    system.F.coeffRef(global_eqn_ids[0], global_var_ids[3 + 2 * i]) = -1.0;
    system.F.coeffRef(global_eqn_ids[i + 1], global_var_ids[3 + 2 * i]) =
        -resistance;
    system.F.coeffRef(global_eqn_ids[i + 1], global_var_ids[0]) = 1.0;
    system.F.coeffRef(global_eqn_ids[i + 1], global_var_ids[2 + 2 * i]) = -1.0;

    system.E.coeffRef(global_eqn_ids[i + 1], global_var_ids[3 + 2 * i]) =
        -inductance;
  }
}

void BloodVesselJunction::update_solution(
    SparseSystem &system, std::vector<double> &parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy) {
  for (size_t i = 0; i < num_outlets; i++) {
    // Get parameters
    auto stenosis_coeff = parameters[global_param_ids[2 * num_outlets + i]];
    auto q_out = y[global_var_ids[3 + 2 * i]];
    auto stenosis_resistance = stenosis_coeff * fabs(q_out);

    // Mass conservation
    system.C(global_eqn_ids[i + 1]) = -stenosis_resistance * q_out;
    system.dC_dy.coeffRef(global_eqn_ids[i + 1], global_var_ids[3 + 2 * i]) =
        -2.0 * stenosis_resistance;
  }
}

void BloodVesselJunction::update_gradient(
    Eigen::SparseMatrix<double> &jacobian,
    Eigen::Matrix<double, Eigen::Dynamic, 1> &residual,
    Eigen::Matrix<double, Eigen::Dynamic, 1> &alpha, std::vector<double> &y,
    std::vector<double> &dy) {
  auto p_in = y[global_var_ids[0]];
  auto q_in = y[global_var_ids[1]];

  residual(global_eqn_ids[0]) = q_in;
  for (size_t i = 0; i < num_outlets; i++) {
    // Get parameters
    auto resistance = alpha[global_param_ids[i]];
    auto inductance = alpha[global_param_ids[num_outlets + i]];
    double stenosis_coeff = 0.0;
    if (global_param_ids.size() / num_outlets > 2) {
      stenosis_coeff = alpha[global_param_ids[2 * num_outlets + i]];
    }
    auto q_out = y[global_var_ids[3 + 2 * i]];
    auto p_out = y[global_var_ids[2 + 2 * i]];
    auto dq_out = dy[global_var_ids[3 + 2 * i]];
    auto stenosis_resistance = stenosis_coeff * fabs(q_out);

    // Resistance
    jacobian.coeffRef(global_eqn_ids[i + 1], global_param_ids[i]) = -q_out;

    // Inductance
    jacobian.coeffRef(global_eqn_ids[i + 1],
                      global_param_ids[num_outlets + i]) = -dq_out;

    // Stenosis Coefficent
    if (global_param_ids.size() / num_outlets > 2) {
      jacobian.coeffRef(global_eqn_ids[i + 1],
                        global_param_ids[2 * num_outlets + i]) =
          -fabs(q_out) * q_out;
    }

    residual(global_eqn_ids[0]) -= q_out;
    residual(global_eqn_ids[i + 1]) =
        p_in - p_out - (resistance + stenosis_resistance) * q_out -
        inductance * dq_out;
  }
}
