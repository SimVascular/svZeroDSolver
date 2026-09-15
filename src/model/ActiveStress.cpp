// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ActiveStress.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

ActiveStress::ActiveStress(
    const std::vector<std::pair<std::string, InputParameter>>& props)
    : input_param_properties(props) {
  for (const auto& p : props) {
    if (p.second.is_number) {
      params_[p.first] = p.second.is_optional ? p.second.default_val : 0.0;
    }
  }
}

void ActiveStress::set_param(const std::string& name, double value) {
  params_[name] = value;
}

std::unique_ptr<ActiveStress> ActiveStress::create(
    const std::string& type_str) {
  if (type_str == "strain_independent") {
    return std::make_unique<StrainIndependentActiveStress>();
  }
  if (type_str == "strain_dependent") {
    return std::make_unique<StrainDependentActiveStress>();
  }
  throw std::runtime_error(
      "Unknown active_stress type '" + type_str +
      "'. Must be one of: strain_independent, strain_dependent");
}

// ============================================================
// StrainIndependentActiveStress
// ============================================================

void StrainIndependentActiveStress::update_constant(
    SparseSystem& system, const std::vector<int>& global_eqn_ids,
    const std::vector<int>& global_var_ids) {
  system.E.coeffRef(global_eqn_ids[NUM_CORE_EQNS], global_var_ids[TAU_VAR]) =
      1;
}

void StrainIndependentActiveStress::update_time(
    SparseSystem& system, const std::vector<int>& global_eqn_ids,
    const std::vector<int>& global_var_ids, double activation_signal) {
  const double alpha_max = params_.at("alpha_max");
  const double alpha_min = params_.at("alpha_min");

  const double act_t =
      alpha_max * activation_signal + alpha_min * (1 - activation_signal);
  act_ = std::abs(act_t);
  act_plus_ = std::max(act_t, 0.0);

  system.F.coeffRef(global_eqn_ids[NUM_CORE_EQNS], global_var_ids[TAU_VAR]) =
      act_;
}

void StrainIndependentActiveStress::update_solution(
    SparseSystem& system, const std::vector<int>& global_eqn_ids,
    const std::vector<int>& global_var_ids,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy, double radius0,
    double activation_signal) {
  (void)global_var_ids;
  (void)y;
  (void)dy;
  (void)radius0;
  (void)activation_signal;

  const double sigma_max = params_.at("sigma_max");
  system.C.coeffRef(global_eqn_ids[NUM_CORE_EQNS]) = -act_plus_ * sigma_max;
}

// ============================================================
// StrainDependentActiveStress
// ============================================================

void StrainDependentActiveStress::update_constant(
    SparseSystem& system, const std::vector<int>& global_eqn_ids,
    const std::vector<int>& global_var_ids) {
  const double alpha_r = params_.at("alpha_r");
  const double mu = params_.at("mu");

  // active stress
  system.F.coeffRef(global_eqn_ids[TAU_EQN], global_var_ids[TAU_VAR]) = -1;

  // relaxation dynamics 
  system.E.coeffRef(global_eqn_ids[OMEGA_EQN], global_var_ids[OMEGA_VAR]) =
      -alpha_r;
  system.F.coeffRef(global_eqn_ids[OMEGA_EQN], global_var_ids[OMEGA_VAR]) =
      -1;

  // sarcomere strain 
  system.E.coeffRef(global_eqn_ids[E_C_EQN], global_var_ids[E_C_VAR]) = -mu;
  system.F.coeffRef(global_eqn_ids[E_C_EQN], global_var_ids[TAU_C_VAR]) = -1;

  // sarcomere stiffness 
  system.E.coeffRef(global_eqn_ids[K_C_EQN], global_var_ids[K_C_VAR]) = -1;

  // sarcomere active stress 
  system.E.coeffRef(global_eqn_ids[TAU_C_EQN], global_var_ids[TAU_C_VAR]) =
      -1;
}

void StrainDependentActiveStress::update_active_stress_values(double e_c,
                                                               double u) {
  // activation strain-dependence
  if (e_c > -0.4 && e_c <= 0.3) {
    n_0_ = 0.22 + 0.53 * e_c;
  } else if (e_c > 0.3 && e_c <= 1.0) {
    n_0_ = 0.112857125 + 0.8871428571 * e_c;
  } else if (e_c > 1.0 && e_c <= 1.3) {
    n_0_ = 1;
  } else if (e_c > 1.3 && e_c <= 2.4) {
    n_0_ = 2.182 - 0.9091 * e_c;
  } else {
    n_0_ = 0.0;
  }

  // relaxation strain-dependence
  if (e_c <= 0.95) {
    m_0_ = 1.87;
  } else if (e_c > 0.95 && e_c <= 1.0) {
    m_0_ = 18.4 - 17.4 * e_c;
  } else {
    m_0_ = 1.0;
  }

  if (u >= 0.0) {
    u_plus_ = u;
    u_minus_ = 0.0;
  } else {
    u_plus_ = 0.0;
    u_minus_ = u;
  }
}

void StrainDependentActiveStress::update_solution(
    SparseSystem& system, const std::vector<int>& global_eqn_ids,
    const std::vector<int>& global_var_ids,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy, double radius0,
    double activation_signal) {
  const double radius = y[global_var_ids[RADIUS_VAR]];
  const double e_c = y[global_var_ids[E_C_VAR]];
  const double tau_c = y[global_var_ids[TAU_C_VAR]];
  const double k_c = y[global_var_ids[K_C_VAR]];
  const double omega = y[global_var_ids[OMEGA_VAR]];
  const double de_c_dt = dy[global_var_ids[E_C_VAR]];

  update_active_stress_values(e_c, activation_signal);

  const double E_s = params_.at("E_s");
  const double alpha = params_.at("alpha");
  const double k_0 = params_.at("k_0");
  const double sigma_0 = params_.at("sigma_0");

  // active stress
  system.C.coeffRef(global_eqn_ids[TAU_EQN]) =
      E_s *
      (-e_c * pow(radius0, 2) + 0.5 * pow(radius, 2) + radius * radius0) /
      (pow(radius0, 2) * pow(2 * e_c + 1, 2));
  system.dC_dy.coeffRef(global_eqn_ids[TAU_EQN], global_var_ids[RADIUS_VAR]) =
      E_s * (1.0 * radius + radius0) / (pow(radius0, 2) * pow(2 * e_c + 1, 2));
  system.dC_dy.coeffRef(global_eqn_ids[TAU_EQN], global_var_ids[E_C_VAR]) =
      E_s *
      (4 * e_c * pow(radius0, 2) - 2.0 * pow(radius, 2) -
       4 * radius * radius0 - pow(radius0, 2) * (2 * e_c + 1)) /
      (pow(radius0, 2) * pow(2 * e_c + 1, 3));

  // relaxation dynamics
  system.C.coeffRef(global_eqn_ids[OMEGA_EQN]) = m_0_;

  // sarcomere element strain
  system.C.coeffRef(global_eqn_ids[E_C_EQN]) =
      E_s * (pow(radius, 2) + 2 * radius * radius0 + pow(radius0, 2)) *
      (-e_c * pow(radius0, 2) + 0.5 * pow(radius, 2) + radius * radius0) /
      (pow(radius0, 4) * pow(2 * e_c + 1, 3));
  system.dC_dy.coeffRef(global_eqn_ids[E_C_EQN], global_var_ids[RADIUS_VAR]) =
      E_s *
      (2 * (radius + radius0) *
           (-e_c * pow(radius0, 2) + 0.5 * pow(radius, 2) +
            radius * radius0) +
       (1.0 * radius + radius0) *
           (pow(radius, 2) + 2 * radius * radius0 + pow(radius0, 2))) /
      (pow(radius0, 4) * pow(2 * e_c + 1, 3));
  system.dC_dy.coeffRef(global_eqn_ids[E_C_EQN], global_var_ids[E_C_VAR]) =
      E_s * (pow(radius, 2) + 2 * radius * radius0 + pow(radius0, 2)) *
      (6 * e_c * pow(radius0, 2) - 3.0 * pow(radius, 2) -
       6 * radius * radius0 - pow(radius0, 2) * (2 * e_c + 1)) /
      (pow(radius0, 4) * pow(2 * e_c + 1, 4));

  // sarcomere stiffness 
  system.C.coeffRef(global_eqn_ids[K_C_EQN]) =
      k_0 * n_0_ * u_plus_ -
      k_c * (alpha * fabs(de_c_dt) + omega * fabs(u_minus_) + u_plus_);
  system.dC_dydot.coeffRef(global_eqn_ids[K_C_EQN], global_var_ids[E_C_VAR]) =
      ((de_c_dt == 0) ? (0) : (-alpha * de_c_dt * k_c / fabs(de_c_dt)));
  system.dC_dy.coeffRef(global_eqn_ids[K_C_EQN], global_var_ids[K_C_VAR]) =
      -alpha * fabs(de_c_dt) - omega * fabs(u_minus_) - u_plus_;
  system.dC_dy.coeffRef(global_eqn_ids[K_C_EQN], global_var_ids[OMEGA_VAR]) =
      -k_c * fabs(u_minus_);

  // sarcomere active stress
  system.C.coeffRef(global_eqn_ids[TAU_C_EQN]) =
      de_c_dt * k_c + n_0_ * sigma_0 * u_plus_ -
      tau_c * (alpha * fabs(de_c_dt) + omega * fabs(u_minus_) + u_plus_);
  system.dC_dydot.coeffRef(global_eqn_ids[TAU_C_EQN],
                           global_var_ids[E_C_VAR]) =
      ((de_c_dt == 0) ? (k_c)
                      : (-alpha * de_c_dt * tau_c / fabs(de_c_dt) + k_c));
  system.dC_dy.coeffRef(global_eqn_ids[TAU_C_EQN], global_var_ids[TAU_C_VAR]) =
      -alpha * fabs(de_c_dt) - omega * fabs(u_minus_) - u_plus_;
  system.dC_dy.coeffRef(global_eqn_ids[TAU_C_EQN], global_var_ids[K_C_VAR]) =
      de_c_dt;
  system.dC_dy.coeffRef(global_eqn_ids[TAU_C_EQN], global_var_ids[OMEGA_VAR]) =
      -tau_c * fabs(u_minus_);
}
