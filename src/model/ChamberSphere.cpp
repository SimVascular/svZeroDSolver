// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ChamberSphere.h"

#include <list>
#include <stdexcept>

#include "Model.h"

namespace {

/// Elastic stress S_el = (1/lambda) Sigma^d : dE/dlambda of the material for
/// the incompressible equibiaxial stretch lambda of the wall,
/// C = diag(lambda^-4, lambda^2, lambda^2) in the (radial, tangential,
/// tangential) basis with the fiber tangential, and its exact derivative from
/// the stress derivative along dC/dlambda = 2 dE/dlambda (J = 1 stays fixed):
/// dS_el/dlambda = -S_el/lambda + (dSigma : E' + Sigma : E'')/lambda with
/// E' = diag(-2 lambda^-5, lambda, lambda), E'' = diag(10 lambda^-6, 1, 1).
double sphere_elastic_stress(const Material& material, double lambda,
                             double& dS_dlambda) {
  using Mat = Material::Mat3;
  const double l2 = lambda * lambda;
  const Mat C = Eigen::Vector3d(1.0 / (l2 * l2), l2, l2).asDiagonal();
  const Mat dE =
      Eigen::Vector3d(-2.0 / (l2 * l2 * lambda), lambda, lambda).asDiagonal();
  const Mat ddE = Eigen::Vector3d(10.0 / (l2 * l2 * l2), 1.0, 1.0).asDiagonal();
  const Material::Vec3 f(0.0, 1.0, 0.0);
  const Mat dC = 2.0 * dE;
  const double dJ = 0.0;
  Mat dS;
  const Mat S = material.stress_derivatives(C, 1.0, f, 1, &dC, &dJ, &dS);
  const double S_el = (S.array() * dE.array()).sum() / lambda;
  dS_dlambda = -S_el / lambda + ((dS.array() * dE.array()).sum() +
                                 (S.array() * ddE.array()).sum()) /
                                    lambda;
  return S_el;
}

}  // namespace

void ChamberSphere::setup_dofs(DOFHandler& dofhandler) {
  if (!activation_function_) {
    throw std::runtime_error("ChamberSphere '" + get_name() +
                             "': activation_function not set. Provide an "
                             "\"activation_function\" JSON block.");
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
  const ActivationSignalType provided_signal =
      activation_function_->signal_type();
  const ActivationSignalType required_signal =
      active_stress_->required_signal_type();
  if (provided_signal != required_signal) {
    throw std::runtime_error(
        "ChamberSphere '" + get_name() + "': active_stress model requires a " +
        (required_signal == ActivationSignalType::rate
             ? "rate-type (e.g. piecewise_rate)"
             : "normalized (e.g. half_cosine)") +
        " activation function, but got a " +
        (provided_signal == ActivationSignalType::rate ? "rate-type"
                                                       : "normalized") +
        " one.");
  }

  std::list<std::string> internal_var_names = {
      "radius", "velo", "stress", "active_stress", "volume", "passive_stress"};
  const std::vector<std::string> states = active_stress_->state_names();
  for (size_t i = 1; i < states.size(); i++) {
    internal_var_names.push_back(states[i]);
  }

  const int num_equations = 7 + active_stress_->num_states();

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

  // passive stress (total stress minus active stress)
  system.F.coeffRef(global_eqn_ids[6], global_var_ids[9]) = -1;
  system.F.coeffRef(global_eqn_ids[6], global_var_ids[6]) = 1;
  system.F.coeffRef(global_eqn_ids[6], global_var_ids[7]) = -1;
}

void ChamberSphere::update_time(SparseSystem& system,
                                std::vector<double>& parameters) {
  (void)system;
  (void)parameters;
  activation_signal_ = activation_function_->compute(model->time);
  active_stress_->update_time(activation_signal_);
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

  // spherical stress: elastic stress of the material plus viscous damping
  const double lambda = 1.0 + radius / radius0;
  double dS_el_dlambda;
  const double S_el = sphere_elastic_stress(*material_, lambda, dS_el_dlambda);
  const double dS_el_dr = dS_el_dlambda / radius0;

  const double C_val_visc = 2 * dradius_dt * eta *
                            (2 * pow(radius0, 12) + pow(radius + radius0, 12)) /
                            (pow(radius0, 2) * pow(radius + radius0, 11));
  const double dC_dy_visc =
      2 * dradius_dt * eta / pow(radius0, 2) -
      44 * dradius_dt * eta * pow(radius0, 10) / pow(radius + radius0, 12);
  const double dC_dydot_visc =
      2 * eta * (2 * pow(radius0, 12) + pow(radius + radius0, 12)) /
      (pow(radius0, 2) * pow(radius + radius0, 11));

  system.C.coeffRef(global_eqn_ids[1]) = S_el + C_val_visc;
  system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[4]) =
      dS_el_dr + dC_dy_visc;
  system.dC_dydot.coeffRef(global_eqn_ids[1], global_var_ids[4]) =
      dC_dydot_visc;

  // volume change
  system.C.coeffRef(global_eqn_ids[2]) =
      4 * M_PI * velo * pow(radius + radius0, 2);
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[4]) =
      8 * M_PI * velo * (radius + radius0);
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[5]) =
      4 * M_PI * pow(radius + radius0, 2);

  // active stress law at the fiber stretch lambda (fibers in the wall plane)
  const int n = active_stress_->num_states();
  std::vector<double> a(n), adot(n);
  for (int i = 0; i < n; i++) {
    a[i] = y[law_var(i)];
    adot[i] = dy[law_var(i)];
  }
  ActiveStress::Linearization lin;
  active_stress_->linearize(a.data(), adot.data(), lambda, lin);
  for (int i = 0; i < n; i++) {
    const int eqn = global_eqn_ids[7 + i];
    system.C.coeffRef(eqn) = lin.R[i];
    for (int j = 0; j < n; j++) {
      system.dC_dy.coeffRef(eqn, law_var(j)) = lin.dR_da[i * n + j];
      system.dC_dydot.coeffRef(eqn, law_var(j)) = lin.dR_dadot[i * n + j];
    }
    system.dC_dy.coeffRef(eqn, global_var_ids[4]) = lin.dR_dlambda[i] / radius0;
  }
}

void ChamberSphere::set_material(std::unique_ptr<Material> m) {
  material_ = std::move(m);
}

void ChamberSphere::set_activation_function(
    std::unique_ptr<ActivationFunction> af) {
  activation_function_ = std::move(af);
}

void ChamberSphere::set_active_stress(std::unique_ptr<ActiveStress> as) {
  active_stress_ = std::move(as);
}
