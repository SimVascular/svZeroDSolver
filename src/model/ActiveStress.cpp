// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ActiveStress.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

// The exact linearizations follow scripts/chamber_tangents.py active_stress
// with scripts/ActiveStress*.yaml.

namespace {

/// Sign of x (0 at x = 0, where |x| has the derivative 0)
inline double sign(double x) { return (x > 0.0) - (x < 0.0); }

}  // namespace

ActiveStress::ActiveStress(
    const std::vector<std::pair<std::string, InputParameter>>& props)
    : ParameterHolder(props) {}

std::unique_ptr<ActiveStress> ActiveStress::create(
    const std::string& type_str) {
  if (type_str == "strain_independent") {
    return std::make_unique<StrainIndependentActiveStress>();
  }
  if (type_str == "strain_dependent") {
    return std::make_unique<StrainDependentActiveStress>();
  }
  if (type_str == "linear_series_spring") {
    return std::make_unique<LinearSeriesSpringActiveStress>();
  }
  throw std::runtime_error(
      "Unknown active_stress type '" + type_str +
      "'. Must be one of: strain_independent, strain_dependent, "
      "linear_series_spring");
}

// ============================================================
// StrainIndependentActiveStress
// ============================================================

void StrainIndependentActiveStress::update_time(double signal) {
  const double alpha_max = params_.at("alpha_max");
  const double alpha_min = params_.at("alpha_min");
  const double act_t = alpha_max * signal + alpha_min * (1 - signal);
  act_ = std::abs(act_t);
  act_plus_ = std::max(act_t, 0.0);
}

void StrainIndependentActiveStress::linearize(const double* a,
                                              const double* adot, double lambda,
                                              Linearization& lin) const {
  (void)lambda;
  lin.reset(1);
  lin.R[0] = adot[0] + act_ * a[0] - act_plus_ * params_.at("sigma_max");
  lin.dR_da[0] = act_;
  lin.dR_dadot[0] = 1.0;
}

// ============================================================
// StrainDependentActiveStress
// ============================================================

void StrainDependentActiveStress::update_time(double signal) {
  if (signal >= 0.0) {
    u_plus_ = signal;
    u_minus_ = 0.0;
  } else {
    u_plus_ = 0.0;
    u_minus_ = signal;
  }
}

void StrainDependentActiveStress::linearize(const double* a, const double* adot,
                                            double lambda,
                                            Linearization& lin) const {
  const double E_s = params_.at("E_s");
  const double mu = params_.at("mu");
  const double alpha_r = params_.at("alpha_r");
  const double alpha = params_.at("alpha");
  const double k_0 = params_.at("k_0");
  const double sigma_0 = params_.at("sigma_0");
  const auto& n0_e = params_vec_.at("n0_e_c");
  const auto& n0_v = params_vec_.at("n0_values");
  const auto& m0_e = params_vec_.at("m0_e_c");
  const auto& m0_v = params_vec_.at("m0_values");
  const double tau = a[0], e_c = a[1], tau_c = a[2], k_c = a[3], omega = a[4],
               sigma_c = a[5], sigma_s = a[6], e_s = a[7];
  const double de_c_dt = adot[1];
  const double n0 = linear_interpolate(e_c, n0_e, n0_v);
  const double dn0 = linear_interpolate_slope(e_c, n0_e, n0_v);
  const double m0 = linear_interpolate(e_c, m0_e, m0_v);
  const double dm0 = linear_interpolate_slope(e_c, m0_e, m0_v);
  const double lam2 = lambda * lambda;
  const double e_1D = 0.5 * (lam2 - 1.0);  // Green-Lagrange fiber strain
  const double g = 1.0 + 2.0 * e_c;
  const double g2 = g * g, g3 = g2 * g;
  const double um = std::abs(u_minus_);
  const double decay = alpha * std::abs(de_c_dt) + omega * um;
  const double sgn = alpha * sign(de_c_dt);
  const int n = 8;
  lin.reset(n);
  // active stress
  lin.R[0] = E_s * (e_1D - e_c) / g2 - tau;
  lin.dR_da[0 * n + 0] = -1.0;
  lin.dR_da[0 * n + 1] = E_s * (g - 2.0 * lam2) / g3;
  lin.dR_dlambda[0] = E_s * lambda / g2;
  // sarcomere element strain
  lin.R[1] =
      E_s * (e_1D - e_c) * (1.0 + 2.0 * e_1D) / g3 - tau_c - mu * de_c_dt;
  lin.dR_da[1 * n + 1] =
      E_s * lam2 * (4.0 * e_c - 3.0 * lam2 + 2.0) / (g2 * g2);
  lin.dR_dadot[1 * n + 1] = -mu;
  lin.dR_da[1 * n + 2] = -1.0;
  lin.dR_dlambda[1] = -E_s * lambda * (g - 2.0 * lam2) / g3;
  // sarcomere active stress
  lin.R[2] = de_c_dt * k_c + n0 * sigma_0 * u_plus_ - tau_c * decay -
             u_plus_ * tau_c - adot[2];
  lin.dR_da[2 * n + 1] = dn0 * sigma_0 * u_plus_;
  lin.dR_dadot[2 * n + 1] = k_c - sgn * tau_c;
  lin.dR_da[2 * n + 2] = -decay - u_plus_;
  lin.dR_dadot[2 * n + 2] = -1.0;
  lin.dR_da[2 * n + 3] = de_c_dt;
  lin.dR_da[2 * n + 4] = -um * tau_c;
  // sarcomere stiffness
  lin.R[3] = k_0 * n0 * u_plus_ - k_c * decay - u_plus_ * k_c - adot[3];
  lin.dR_da[3 * n + 1] = dn0 * k_0 * u_plus_;
  lin.dR_dadot[3 * n + 1] = -sgn * k_c;
  lin.dR_da[3 * n + 3] = -decay - u_plus_;
  lin.dR_dadot[3 * n + 3] = -1.0;
  lin.dR_da[3 * n + 4] = -um * k_c;
  // relaxation dynamics
  lin.R[4] = m0 - omega - alpha_r * adot[4];
  lin.dR_da[4 * n + 1] = dm0;
  lin.dR_da[4 * n + 4] = -1.0;
  lin.dR_dadot[4 * n + 4] = -alpha_r;
  // sarcomere stress
  lin.R[5] = mu * de_c_dt + tau_c - sigma_c;
  lin.dR_dadot[5 * n + 1] = mu;
  lin.dR_da[5 * n + 2] = 1.0;
  lin.dR_da[5 * n + 5] = -1.0;
  // series element stress
  lin.R[6] = tau * g - sigma_s;
  lin.dR_da[6 * n + 0] = g;
  lin.dR_da[6 * n + 1] = 2.0 * tau;
  lin.dR_da[6 * n + 6] = -1.0;
  // series element strain
  lin.R[7] = sigma_s / E_s - e_s;
  lin.dR_da[7 * n + 6] = 1.0 / E_s;
  lin.dR_da[7 * n + 7] = -1.0;
}

// ============================================================
// LinearSeriesSpringActiveStress
// ============================================================

void LinearSeriesSpringActiveStress::finalize() {
  // Default Frank-Starling curve n0(e_c): that of the genet23 authors' input
  // functions (genet23 Table 1 cites it from its Ref. 6)
  if (params_vec_.at("n0_e_c").size() < 2 &&
      params_vec_.at("n0_values").size() < 2) {
    params_vec_["n0_e_c"] = {-0.1668, -0.0073, 0.0534, 0.0969, 0.1326,
                             0.2016,  0.4663,  0.9187, 1.1762};
    params_vec_["n0_values"] = {0.0, 0.5614, 0.7748, 0.8933, 0.9618,
                                1.0, 1.0,    0.1075, 0.0};
  }
  if (params_vec_.at("n0_e_c").size() != params_vec_.at("n0_values").size() ||
      params_vec_.at("n0_e_c").size() < 2) {
    throw std::runtime_error(
        "linear_series_spring active_stress: n0_e_c and n0_values must have "
        "equal length "
        ">= 2.");
  }
  if (params_.at("mu") <= 0.0) {
    throw std::runtime_error(
        "linear_series_spring active_stress: mu must be positive.");
  }
}

void LinearSeriesSpringActiveStress::update_time(double signal) {
  nu_abs_ = std::abs(signal);
  nu_plus_ = std::max(signal, 0.0);
}

void LinearSeriesSpringActiveStress::linearize(const double* a,
                                               const double* adot,
                                               double lambda,
                                               Linearization& lin) const {
  const double k_s = params_.at("k_s");
  const double k_0 = params_.at("k_0");
  const double mu = params_.at("mu");
  const double alpha = params_.at("alpha");
  const double sigma_0 = params_.at("sigma_0");
  const auto& n0_e = params_vec_.at("n0_e_c");
  const auto& n0_v = params_vec_.at("n0_values");
  const double tau = a[0], e_c = a[1], tau_c = a[2], k_c = a[3];
  const double ecdot = adot[1];
  const double n0 = linear_interpolate(e_c, n0_e, n0_v);
  const double dn0 = linear_interpolate_slope(e_c, n0_e, n0_v);
  const double decay = nu_abs_ + alpha * std::abs(ecdot);
  const double sgn = alpha * sign(ecdot);
  const double s1d = tau_c + mu * ecdot;  // 1D fiber stress
  const int n = 4;
  lin.reset(n);
  // active stress tau = (tau_c + mu e_c_dot) / (1 + e_fib)
  lin.R[0] = tau - s1d / lambda;
  lin.dR_da[0 * n + 0] = 1.0;
  lin.dR_dadot[0 * n + 1] = -mu / lambda;
  lin.dR_da[0 * n + 2] = -1.0 / lambda;
  lin.dR_dlambda[0] = s1d / (lambda * lambda);
  // series-spring force balance: mu * e_c_dot = k_s (e_fib - e_c) - tau_c
  lin.R[1] = ecdot - (k_s * (lambda - 1.0 - e_c) - tau_c) / mu;
  lin.dR_da[1 * n + 1] = k_s / mu;
  lin.dR_dadot[1 * n + 1] = 1.0;
  lin.dR_da[1 * n + 2] = 1.0 / mu;
  lin.dR_dlambda[1] = -k_s / mu;
  // tau_c dot = -decay*tau_c + e_c_dot*k_c + n0*sigma0*nu_+
  lin.R[2] = adot[2] - (-decay * tau_c + ecdot * k_c + n0 * sigma_0 * nu_plus_);
  lin.dR_da[2 * n + 1] = -dn0 * sigma_0 * nu_plus_;
  lin.dR_dadot[2 * n + 1] = sgn * tau_c - k_c;
  lin.dR_da[2 * n + 2] = decay;
  lin.dR_dadot[2 * n + 2] = 1.0;
  lin.dR_da[2 * n + 3] = -ecdot;
  // k_c dot = -decay*k_c + n0*k_0*nu_+
  lin.R[3] = adot[3] - (-decay * k_c + n0 * k_0 * nu_plus_);
  lin.dR_da[3 * n + 1] = -dn0 * k_0 * nu_plus_;
  lin.dR_dadot[3 * n + 1] = sgn * k_c;
  lin.dR_da[3 * n + 3] = decay;
  lin.dR_dadot[3 * n + 3] = 1.0;
}
