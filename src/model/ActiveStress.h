// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

/**
 * @file ActiveStress.h
 * @brief Active stress laws shared by the chamber blocks
 */

#ifndef SVZERODSOLVER_MODEL_ACTIVESTRESS_HPP_
#define SVZERODSOLVER_MODEL_ACTIVESTRESS_HPP_

#include <complex>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "ActivationFunction.h"
#include "Parameter.h"
#include "ParameterHolder.h"

/**
 * @brief Active stress law at a material point of a chamber wall
 *
 * A law defines the active fiber stress \f$\tau\f$, the 2nd Piola-Kirchhoff
 * stress \f$\Sigma^a = \tau\,\mathbf{e}_F\otimes\mathbf{e}_F\f$ along the
 * fiber direction, through a set of states \f$a\f$ (state 0 is \f$\tau\f$
 * itself) and as many equations
 * \f[
 * R(a, \dot a, \lambda) = 0,
 * \f]
 * where \f$\lambda = \sqrt{\mathbf{e}_F\cdot C\cdot\mathbf{e}_F}\f$ is the
 * fiber stretch. The activation signal of the chamber's \ref
 * ActivationFunction is set once per time step with update_time(). The law is the same at every material point,
 * so a block evaluates it once for a uniform wall (\ref ChamberSphere) or at
 * each integration point (\ref ChamberCylinder).
 *
 * Types (JSON `active_stress` object): `strain_independent`,
 * `strain_dependent`, `genet23`.
 */
class ActiveStress : public ParameterHolder {
 public:
  using cplx = std::complex<double>;  ///< complex-step scalar

  /**
   * @brief Construct an ActiveStress law
   *
   * @param props Parameter name/spec pairs for this active stress type
   */
  explicit ActiveStress(
      const std::vector<std::pair<std::string, InputParameter>>& props);

  /**
   * @brief Number of states (and equations) at a material point
   *
   * @return Number of states, state 0 being the active fiber stress
   */
  virtual int num_states() const = 0;

  /**
   * @brief Names of the states, state 0 being `active_stress`
   *
   * @return State names in order
   */
  virtual std::vector<std::string> state_names() const = 0;

  /**
   * @brief Type of activation signal this law needs
   *
   * Checked by the chamber blocks against ActivationFunction::signal_type.
   *
   * @return The required ActivationSignalType
   */
  virtual ActivationSignalType required_signal_type() const {
    return ActivationSignalType::normalized;
  }

  /**
   * @brief Set the activation signal of the current time
   *
   * @param signal Activation signal of the chamber's ActivationFunction
   */
  virtual void update_time(double signal) = 0;

  /**
   * @brief Residuals of the state equations at a material point
   *
   * @param a States (num_states())
   * @param adot State rates
   * @param lambda Fiber stretch
   * @param R Residuals (num_states())
   */
  virtual void residual(const double* a, const double* adot, double lambda,
                        double* R) const = 0;

  /**
   * @brief Residuals for complex arguments (complex-step tangents)
   *
   * @param a States
   * @param adot State rates
   * @param lambda Fiber stretch
   * @param R Residuals
   */
  virtual void residual(const cplx* a, const cplx* adot, cplx lambda,
                        cplx* R) const = 0;

  /**
   * @brief Residuals and their derivatives at a material point
   */
  struct Linearization {
    std::vector<double> R;           ///< residuals
    std::vector<double> dR_da;       ///< dR_i/da_j at [i * n + j]
    std::vector<double> dR_dadot;    ///< dR_i/dadot_j at [i * n + j]
    std::vector<double> dR_dlambda;  ///< dR_i/dlambda
  };

  /**
   * @brief Residuals and their derivatives by complex step
   *
   * @param a States
   * @param adot State rates
   * @param lambda Fiber stretch
   * @param lin Residuals and derivatives
   */
  void linearize(const double* a, const double* adot, double lambda,
                 Linearization& lin) const;

  /**
   * @brief Factory: create an active stress law from a type string
   *
   * @param type_str One of: "strain_independent", "strain_dependent",
   * "genet23"
   * @return Unique pointer to the created active stress law
   */
  static std::unique_ptr<ActiveStress> create(const std::string& type_str);
};

/**
 * @brief Strain-independent active stress (\cite pfaller2019importance)
 *
 * One state, the active stress \f$\tau\f$:
 * \f[
 * \dot{\tau} + |a| \tau - \sigma_\text{max} a_+ = 0, \quad a_+ = \max(a, 0),
 * \quad a = f \alpha_\text{max} + (1 - f) \alpha_\text{min}
 * \f]
 * where \f$f \in [0, 1]\f$ is the (normalized) activation function.
 *
 * Parameters: `alpha_max`, `alpha_min`, `sigma_max`
 */
class StrainIndependentActiveStress : public ActiveStress {
 public:
  StrainIndependentActiveStress()
      : ActiveStress({{"alpha_max", InputParameter()},
                      {"alpha_min", InputParameter()},
                      {"sigma_max", InputParameter()}}) {}

  int num_states() const override { return 1; }
  std::vector<std::string> state_names() const override {
    return {"active_stress"};
  }
  void update_time(double signal) override;
  void residual(const double* a, const double* adot, double lambda,
                double* R) const override;
  void residual(const cplx* a, const cplx* adot, cplx lambda,
                cplx* R) const override;

 private:
  double act_ = 0.0;       ///< |a|
  double act_plus_ = 0.0;  ///< a_+
};

/**
 * @brief Strain-dependent active stress model from \cite caruel13
 *
 * States \f$\tau, e_c, \tau_c, k_c, \omega, \sigma_c, \sigma_s, e_s\f$:
 * \f[
 * \tau = E_s \frac{e_{1D} - e_c}{(1 + 2 e_c)^2}
 * \f]
 * with the Green-Lagrange fiber strain \f$e_{1D} = (\lambda^2 - 1)/2\f$, and
 * evolution equations
 * \f[
 * \dot{\omega} = \frac{1}{\alpha_r}(m_0 - \omega)
 * \f]
 * \f[
 * \dot{e}_c = \frac{1}{\mu}\left(E_s \frac{(e_{1D} - e_c)(1 + 2 e_{1D})}
 * {(1 + 2 e_c)^3} - \tau_c\right)
 * \f]
 * \f[
 * \dot{k}_c = -(|u|_+ + \omega |u|_- + \alpha |\dot{e}_c|) k_c + n_0 k_0 |u|_+
 * \f]
 * \f[
 * \dot{\tau}_c = -(|u|_+ + \omega |u|_- + \alpha |\dot{e}_c|) \tau_c + n_0
 * \sigma_0 |u|_+ + k_c \dot{e}_c
 * \f]
 * and the additional outputs (sarcomere stress, series element stress and
 * strain)
 * \f[
 * \sigma_c = \tau_c + \mu \dot{e}_c, \quad
 * \sigma_s = \tau (1 + 2 e_c), \quad
 * e_s = \frac{\sigma_s}{E_s}
 * \f]
 * where \f$n_0\f$ (activation strain-dependence) and \f$m_0\f$ (relaxation
 * strain-dependence) are piecewise-linear functions of \f$e_c\f$, given as
 * break point tables (`n0_e_c`/`n0_values` and `m0_e_c`/`m0_values`) and
 * evaluated with \ref linear_interpolate, and \f$u(t)\f$ is a reaction-rate
 * signal supplied by the chamber's ActivationFunction (e.g. \ref
 * PiecewiseRateActivation).
 *
 * Parameters: `E_s`, `mu`, `alpha_r`, `alpha`, `k_0`, `sigma_0`, `n0_e_c`,
 * `n0_values`, `m0_e_c`, `m0_values`
 */
class StrainDependentActiveStress : public ActiveStress {
 public:
  StrainDependentActiveStress()
      : ActiveStress({{"E_s", InputParameter()},
                      {"mu", InputParameter()},
                      {"alpha_r", InputParameter()},
                      {"alpha", InputParameter()},
                      {"k_0", InputParameter()},
                      {"sigma_0", InputParameter()},
                      {"n0_e_c", InputParameter(false, true)},
                      {"n0_values", InputParameter(false, true)},
                      {"m0_e_c", InputParameter(false, true)},
                      {"m0_values", InputParameter(false, true)}}) {}

  int num_states() const override { return 8; }
  std::vector<std::string> state_names() const override {
    return {"active_stress", "e_c",     "tau_c",   "k_c",
            "omega",         "sigma_c", "sigma_s", "e_s"};
  }
  ActivationSignalType required_signal_type() const override {
    return ActivationSignalType::rate;
  }
  void update_time(double signal) override;
  void residual(const double* a, const double* adot, double lambda,
                double* R) const override;
  void residual(const cplx* a, const cplx* adot, cplx lambda,
                cplx* R) const override;

 private:
  template <typename T>
  void residual_t(const T* a, const T* adot, T lambda, T* R) const;

  double u_plus_ = 0.0;   ///< positive part of the reaction-rate signal
  double u_minus_ = 0.0;  ///< negative part of the reaction-rate signal
};

/**
 * @brief Bestel-Clement-Sorine active stress of \cite genet23
 *
 * The law of \cite genet23 (Eqs. 30-33, after Chapelle et al. 2012) with the
 * fiber strain \f$e_\text{fib} = \lambda - 1\f$ and a linear series spring.
 * States \f$\tau, e_c, \tau_c, k_c\f$:
 * \f[
 * \tau = \frac{\tau_c + \mu\dot e_c}{\lambda},\quad
 * \mu\dot e_c = k_s(e_\text{fib}-e_c)-\tau_c,
 * \f]
 * \f[
 * \dot\tau_c = -(|\nu|+\alpha|\dot e_c|)\tau_c + \dot e_c\,k_c +
 *   n_0(e_c)\sigma_0|\nu|_+,\quad
 * \dot k_c = -(|\nu|+\alpha|\dot e_c|)k_c + n_0(e_c)k_0|\nu|_+,
 * \f]
 * driven by a rate-type activation signal \f$\nu(t)\f$ (e.g. the ECG-derived
 * activation rate of the authors' input functions as a `piecewise_rate`
 * table). Writing \f$\tau\f$ through
 * \f$\tau_c+\mu\dot e_c\f$ instead of the equivalent
 * \f$k_s(e_\text{fib}-e_c)\f$ keeps the stiff \f$k_s\f$ out of the wall
 * stress. The Frank-Starling factor \f$n_0(e_c)\f$ is a piecewise-linear
 * table (`n0_e_c`/`n0_values`), by default that of the genet23 authors' code.
 *
 * Parameters: `sigma_0`, `k_s` (default 1e8), `k_0` (default 260e3), `mu`
 * (default 70), `alpha` (default 0), `n0_e_c`, `n0_values` (optional)
 */
class Genet23ActiveStress : public ActiveStress {
 public:
  Genet23ActiveStress()
      : ActiveStress({{"sigma_0", InputParameter()},
                      {"k_s", InputParameter(true, false, true, 1.0e8)},
                      {"k_0", InputParameter(true, false, true, 260.0e3)},
                      {"mu", InputParameter(true, false, true, 70.0)},
                      {"alpha", InputParameter(true, false, true, 0.0)},
                      {"n0_e_c", InputParameter(true, true, true, 0.0)},
                      {"n0_values", InputParameter(true, true, true, 0.0)}}) {}

  void finalize() override;
  int num_states() const override { return 4; }
  std::vector<std::string> state_names() const override {
    return {"active_stress", "e_c", "tau_c", "k_c"};
  }
  ActivationSignalType required_signal_type() const override {
    return ActivationSignalType::rate;
  }
  void update_time(double signal) override;
  void residual(const double* a, const double* adot, double lambda,
                double* R) const override;
  void residual(const cplx* a, const cplx* adot, cplx lambda,
                cplx* R) const override;

 private:
  template <typename T>
  void residual_t(const T* a, const T* adot, T lambda, T* R) const;

  double nu_abs_ = 0.0;   ///< |nu(t)|
  double nu_plus_ = 0.0;  ///< |nu(t)|_+
};

#endif  // SVZERODSOLVER_MODEL_ACTIVESTRESS_HPP_
