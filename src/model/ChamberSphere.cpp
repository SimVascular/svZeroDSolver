// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ChamberSphere.h"

#include <list>
#include <stdexcept>

#include "Model.h"

void ChamberSphere::setup_dofs(DOFHandler& dofhandler) {
  if(!activation_function_) {
    throw std::runtime_error(
        "ChamberSphere '" + get_name() +
        "': activation_function not set. Provide an \"activation_function\" JSON block.");
  }
  if (!active_stress_) {
    throw std::runtime_error(
        "ChamberSphere '" + get_name() +
        "': active_stress not set. Provide an \"active_stress\" JSON block.");
  }
  if (!material_) {
    throw std::runtime_error(
        "ChamberSphere '" + get_name() +
        "': material not set. Provide a \"material\" JSON block.");
  }

  std::list<std::string> internal_var_names = {"radius", "velo", "stress",
                                               "tau", "volume"};
  for (const auto& name : active_stress_->extra_var_names()) {
    internal_var_names.push_back(name);
  }

  const int num_equations =
      ActiveStress::NUM_CORE_EQNS + active_stress_->num_extra_equations();

  Block::setup_dofs_(dofhandler, num_equations, internal_var_names);
}

void ChamberSphere::update_constant(SparseSystem& system,
                                    std::vector<double>& parameters) {
  const double rho = parameters[global_param_ids[ParamId::rho]];
  const double thick0 = parameters[global_param_ids[ParamId::thick0]];

  // balance of linear momentum
  system.E.coeffRef(global_eqn_ids[0], global_var_ids[5]) = rho * thick0;

  // spherical stress
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[6]) = -1;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[7]) = 1;

  // volume change
  system.E.coeffRef(global_eqn_ids[2], global_var_ids[8]) = -1;

  // acceleration
  system.E.coeffRef(global_eqn_ids[3], global_var_ids[4]) = 1;
  system.F.coeffRef(global_eqn_ids[3], global_var_ids[5]) = -1;

  // conservation of mass
  system.F.coeffRef(global_eqn_ids[4], global_var_ids[1]) = 1;
  system.F.coeffRef(global_eqn_ids[4], global_var_ids[3]) = -1;
  system.E.coeffRef(global_eqn_ids[4], global_var_ids[8]) = -1;

  // pressure equality
  system.F.coeffRef(global_eqn_ids[5], global_var_ids[0]) = 1;
  system.F.coeffRef(global_eqn_ids[5], global_var_ids[2]) = -1;

  // active stress
  active_stress_->update_constant(system, global_eqn_ids, global_var_ids);
}

void ChamberSphere::update_time(SparseSystem& system,
                                std::vector<double>& parameters) {
  activation_signal_ = activation_function_->compute(model->time);
  active_stress_->update_time(system, global_eqn_ids, global_var_ids,
                              activation_signal_);
}

void ChamberSphere::update_solution(
    SparseSystem& system, std::vector<double>& parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) {
  const double thick0 = parameters[global_param_ids[ParamId::thick0]];
  const double radius0 = parameters[global_param_ids[ParamId::radius0]];
  const double eta = parameters[global_param_ids[ParamId::eta]];

  const double velo = y[global_var_ids[5]];
  const double Pout = y[global_var_ids[2]];
  const double stress = y[global_var_ids[6]];
  const double dradius_dt = dy[global_var_ids[4]];
  const double radius = y[global_var_ids[4]];

  // balance of momentum
  system.C.coeffRef(global_eqn_ids[0]) =
      (radius + radius0) * (-Pout * (radius + radius0) + stress * thick0) /
      pow(radius0, 2);
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[2]) =
      -pow(radius + radius0, 2) / pow(radius0, 2);
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[4]) =
      (-2 * Pout * (radius + radius0) + stress * thick0) / pow(radius0, 2);
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[6]) =
      thick0 * (radius + radius0) / pow(radius0, 2);

  // spherical stress: material-dependent elastic term plus viscous damping
  const auto mat = material_->compute(radius, radius0);

  const double C_val_visc =
      2 * dradius_dt * eta * (2 * pow(radius0, 12) + pow(radius + radius0, 12)) /
      (pow(radius0, 2) * pow(radius + radius0, 11));
  const double dC_dy_visc =
      2 * dradius_dt * eta / pow(radius0, 2) -
      44 * dradius_dt * eta * pow(radius0, 10) / pow(radius + radius0, 12);
  const double dC_dydot_visc =
      2 * eta * (2 * pow(radius0, 12) + pow(radius + radius0, 12)) /
      (pow(radius0, 2) * pow(radius + radius0, 11));

  system.C.coeffRef(global_eqn_ids[1]) = mat.C_val + C_val_visc;
  system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[4]) =
      mat.dC_dy_radius + dC_dy_visc;
  system.dC_dydot.coeffRef(global_eqn_ids[1], global_var_ids[4]) =
      dC_dydot_visc;

  // volume change
  system.C.coeffRef(global_eqn_ids[2]) =
      4 * M_PI * velo * pow(radius + radius0, 2);
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[4]) =
      8 * M_PI * velo * (radius + radius0);
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[5]) =
      4 * M_PI * pow(radius + radius0, 2);

  // active stress
  active_stress_->update_solution(system, global_eqn_ids, global_var_ids, y,
                                  dy, radius0, activation_signal_);
}

void ChamberSphere::set_material(std::unique_ptr<SphereMaterial> m) {
  material_ = std::move(m);
}

void ChamberSphere::set_activation_function(
    std::unique_ptr<ActivationFunction> af) {
  activation_function_ = std::move(af);
}

void ChamberSphere::set_active_stress(std::unique_ptr<ActiveStress> as) {
  active_stress_ = std::move(as);
}
