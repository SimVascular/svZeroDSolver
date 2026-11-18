// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ActiveStress.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {

using cplx = std::complex<double>;

/// Real part of a scalar, for branches on complex-step arguments.
inline double re(double x) { return x; }
inline double re(const cplx& x) { return x.real(); }

/// |x| compatible with complex-step, with derivative 0 at x = 0.
template <typename T>
inline T cs_abs(const T& x) {
  if (re(x) > 0.0) return x;
  if (re(x) < 0.0) return -x;
  return T(0.0);
}

/// Piecewise-linear table y(x) (\ref linear_interpolate) compatible with
/// complex-step: the value at re(x) plus the slope times the perturbation.
template <typename T>
inline T table(const T& x, const std::vector<double>& xs,
               const std::vector<double>& ys) {
  const double x0 = re(x);
  return linear_interpolate(x0, xs, ys) +
         linear_interpolate_slope(x0, xs, ys) * (x - x0);
}

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

void ActiveStress::linearize_complex_step(const double* a, const double* adot,
                                          double lambda,
                                          Linearization& lin) const {
  const int n = num_states();
  lin.R.assign(n, 0.0);
  lin.dR_da.assign(n * n, 0.0);
  lin.dR_dadot.assign(n * n, 0.0);
  lin.dR_dlambda.assign(n, 0.0);
  residual(a, adot, lambda, lin.R.data());

  const double h = 1e-30;
  std::vector<cplx> ac(a, a + n), adc(adot, adot + n), Rc(n);
  for (int j = 0; j < n; j++) {
    ac[j] += cplx(0.0, h);
    residual(ac.data(), adc.data(), cplx(lambda), Rc.data());
    for (int i = 0; i < n; i++) lin.dR_da[i * n + j] = std::imag(Rc[i]) / h;
    ac[j] = a[j];
    adc[j] += cplx(0.0, h);
    residual(ac.data(), adc.data(), cplx(lambda), Rc.data());
    for (int i = 0; i < n; i++) lin.dR_dadot[i * n + j] = std::imag(Rc[i]) / h;
    adc[j] = adot[j];
  }
  residual(ac.data(), adc.data(), cplx(lambda, h), Rc.data());
  for (int i = 0; i < n; i++) lin.dR_dlambda[i] = std::imag(Rc[i]) / h;
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

void StrainIndependentActiveStress::residual(const double* a,
                                             const double* adot, double lambda,
                                             double* R) const {
  (void)lambda;
  R[0] = adot[0] + act_ * a[0] - act_plus_ * params_.at("sigma_max");
}

void StrainIndependentActiveStress::residual(const cplx* a, const cplx* adot,
                                             cplx lambda, cplx* R) const {
  (void)lambda;
  R[0] = adot[0] + act_ * a[0] - act_plus_ * params_.at("sigma_max");
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

template <typename T>
void StrainDependentActiveStress::residual_t(const T* a, const T* adot,
                                             T lambda, T* R) const {
  const double E_s = params_.at("E_s");
  const double mu = params_.at("mu");
  const double alpha_r = params_.at("alpha_r");
  const double alpha = params_.at("alpha");
  const double k_0 = params_.at("k_0");
  const double sigma_0 = params_.at("sigma_0");
  const T &tau = a[0], &e_c = a[1], &tau_c = a[2], &k_c = a[3], &omega = a[4],
          &sigma_c = a[5], &sigma_s = a[6], &e_s = a[7];
  const T& de_c_dt = adot[1];

  const T n_0 =
      table(e_c, params_vec_.at("n0_e_c"), params_vec_.at("n0_values"));
  const T m_0 =
      table(e_c, params_vec_.at("m0_e_c"), params_vec_.at("m0_values"));
  const T e_1D = 0.5 * (lambda * lambda - 1.0);  // Green-Lagrange fiber strain
  const T one_2ec = 1.0 + 2.0 * e_c;
  const T decay = alpha * cs_abs(de_c_dt) + omega * std::abs(u_minus_);

  // active stress
  R[0] = E_s * (e_1D - e_c) / (one_2ec * one_2ec) - tau;
  // sarcomere element strain
  R[1] =
      E_s * (e_1D - e_c) * (1.0 + 2.0 * e_1D) / (one_2ec * one_2ec * one_2ec) -
      tau_c - mu * de_c_dt;
  // sarcomere active stress
  R[2] = de_c_dt * k_c + n_0 * sigma_0 * u_plus_ - tau_c * decay -
         u_plus_ * tau_c - adot[2];
  // sarcomere stiffness
  R[3] = k_0 * n_0 * u_plus_ - k_c * decay - u_plus_ * k_c - adot[3];
  // relaxation dynamics
  R[4] = m_0 - omega - alpha_r * adot[4];
  // sarcomere stress
  R[5] = mu * de_c_dt + tau_c - sigma_c;
  // series element stress
  R[6] = tau * one_2ec - sigma_s;
  // series element strain
  R[7] = sigma_s / E_s - e_s;
}

void StrainDependentActiveStress::residual(const double* a, const double* adot,
                                           double lambda, double* R) const {
  residual_t<double>(a, adot, lambda, R);
}

void StrainDependentActiveStress::residual(const cplx* a, const cplx* adot,
                                           cplx lambda, cplx* R) const {
  residual_t<cplx>(a, adot, lambda, R);
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

template <typename T>
void LinearSeriesSpringActiveStress::residual_t(const T* a, const T* adot,
                                                T lambda, T* R) const {
  const double k_s = params_.at("k_s");
  const double k_0 = params_.at("k_0");
  const double mu = params_.at("mu");
  const double alpha = params_.at("alpha");
  const double sigma_0 = params_.at("sigma_0");
  const T &e_c = a[1], &tau_c = a[2], &k_c = a[3];
  const T& ecdot = adot[1];
  const T n0 =
      table(e_c, params_vec_.at("n0_e_c"), params_vec_.at("n0_values"));
  const T decay = nu_abs_ + alpha * cs_abs(ecdot);
  const T e_fib = lambda - 1.0;
  // active stress tau = (tau_c + mu e_c_dot) / (1 + e_fib)
  R[0] = a[0] - (tau_c + mu * ecdot) / lambda;
  // series-spring force balance: mu * e_c_dot = k_s (e_fib - e_c) - tau_c
  R[1] = ecdot - (k_s * (e_fib - e_c) - tau_c) / mu;
  // tau_c dot = -decay*tau_c + e_c_dot*k_c + n0*sigma0*nu_+
  R[2] = adot[2] - (-decay * tau_c + ecdot * k_c + n0 * sigma_0 * nu_plus_);
  // k_c dot = -decay*k_c + n0*k_0*nu_+
  R[3] = adot[3] - (-decay * k_c + n0 * k_0 * nu_plus_);
}

void LinearSeriesSpringActiveStress::residual(const double* a,
                                              const double* adot, double lambda,
                                              double* R) const {
  residual_t<double>(a, adot, lambda, R);
}

void LinearSeriesSpringActiveStress::residual(const cplx* a, const cplx* adot,
                                              cplx lambda, cplx* R) const {
  residual_t<cplx>(a, adot, lambda, R);
}

// ============================================================
// Exact linearizations (scripts/chamber_tangents.py active_stress with
// scripts/ActiveStress*.yaml)
// ============================================================

namespace {

/// Sign of x (0 at x = 0, where |x| has the derivative 0)
inline double sign(double x) { return (x > 0.0) - (x < 0.0); }

}  // namespace

void StrainIndependentActiveStress::linearize(const double* a,
                                              const double* adot, double lambda,
                                              Linearization& lin) const {
  (void)lambda;
  lin.reset(1);
  lin.R[0] = adot[0] + act_ * a[0] - act_plus_ * params_.at("sigma_max");
  lin.dR_da[0] = act_;
  lin.dR_dadot[0] = 1.0;
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
