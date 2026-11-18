// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

/**
 * @file ActiveStress.h
 * @brief Active stress laws shared by the chamber blocks
 */

#ifndef SVZERODSOLVER_MODEL_ACTIVESTRESS_HPP_
#define SVZERODSOLVER_MODEL_ACTIVESTRESS_HPP_

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
 * ActivationFunction is set once per time step with update_time(). The law is
 * the same at every material point, so a block evaluates it once for a uniform
 * wall (\ref ChamberSphere) or at each integration point (\ref
 * ChamberCylinder).
 *
 * ### Linearization
 *
 * linearize() returns the residuals \f$R\f$ with their exact derivatives
 * \f$\partial R/\partial a\f$, \f$\partial R/\partial\dot a\f$ and
 * \f$\partial R/\partial\lambda\f$. The blocks chain the last with the
 * derivative of the fiber stretch,
 * \f$\partial\lambda/\partial C =
 * \mathbf{e}_F\otimes\mathbf{e}_F/(2\lambda)\f$. The derivatives follow from
 * the residuals in `scripts/ActiveStress*.yaml` with
 * `scripts/chamber_tangents.py active_stress`. A piecewise-linear table, such
 * as \f$n_0(e_c)\f$, contributes its local slope \f$n_0'(e_c)\f$, and
 * \f$|x|\f$ the derivative \f$\operatorname{sign}(x)\f$ (0 at
 * \f$x = 0\f$). Each law below lists its residuals in the order of its
 * states and their nonzero derivatives.
 *
 * Types (JSON `active_stress` object): `strain_independent`,
 * `strain_dependent`, `linear_series_spring`.
 */
class ActiveStress : public ParameterHolder {
 public:
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
   * @brief Residuals and their derivatives at a material point
   */
  struct Linearization {
    std::vector<double> R;           ///< residuals
    std::vector<double> dR_da;       ///< dR_i/da_j at [i * n + j]
    std::vector<double> dR_dadot;    ///< dR_i/dadot_j at [i * n + j]
    std::vector<double> dR_dlambda;  ///< dR_i/dlambda

    /**
     * @brief Zero all entries for n states
     *
     * @param n Number of states
     */
    void reset(int n) {
      R.assign(n, 0.0);
      dR_da.assign(n * n, 0.0);
      dR_dadot.assign(n * n, 0.0);
      dR_dlambda.assign(n, 0.0);
    }
  };

  /**
   * @brief Residuals and their exact derivatives at a material point
   *
   * @param a States
   * @param adot State rates
   * @param lambda Fiber stretch
   * @param lin Residuals and derivatives
   */
  virtual void linearize(const double* a, const double* adot, double lambda,
                         Linearization& lin) const = 0;

  /**
   * @brief Factory: create an active stress law from a type string
   *
   * @param type_str One of: "strain_independent", "strain_dependent",
   * "linear_series_spring"
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
 * Linearization of \f$R = \dot\tau + |a|\tau - \sigma_\text{max} a_+\f$:
 * \f$\partial R/\partial\tau = |a|\f$,
 * \f$\partial R/\partial\dot\tau = 1\f$; it does not depend on
 * \f$\lambda\f$.
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
  void linearize(const double* a, const double* adot, double lambda,
                 Linearization& lin) const override;

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
 * Linearization, with \f$g = 1 + 2e_c\f$ and the decay
 * \f$d = \alpha|\dot e_c| + \omega|u|_-\f$:
 * - \f$R_0 = E_s(e_{1D} - e_c)/g^2 - \tau\f$:
 *   \f$\partial_\tau R_0 = -1\f$,
 *   \f$\partial_{e_c} R_0 = E_s(g - 2\lambda^2)/g^3\f$,
 *   \f$\partial_\lambda R_0 = E_s\lambda/g^2\f$;
 * - \f$R_1 = E_s(e_{1D} - e_c)(1 + 2e_{1D})/g^3 - \tau_c - \mu\dot e_c\f$:
 *   \f$\partial_{e_c} R_1 = E_s\lambda^2(4e_c - 3\lambda^2 + 2)/g^4\f$,
 *   \f$\partial_{\tau_c} R_1 = -1\f$, \f$\partial_{\dot e_c} R_1 = -\mu\f$,
 *   \f$\partial_\lambda R_1 = -E_s\lambda(g - 2\lambda^2)/g^3\f$;
 * - \f$R_2 = \dot e_c k_c + n_0\sigma_0|u|_+ - (d + |u|_+)\tau_c -
 * \dot\tau_c\f$:
 *   \f$\partial_{e_c} R_2 = n_0'\sigma_0|u|_+\f$,
 *   \f$\partial_{\dot e_c} R_2 = k_c - \alpha\operatorname{sign}(\dot
 * e_c)\tau_c\f$,
 *   \f$\partial_{\tau_c} R_2 = -d - |u|_+\f$,
 *   \f$\partial_{\dot\tau_c} R_2 = -1\f$,
 *   \f$\partial_{k_c} R_2 = \dot e_c\f$,
 *   \f$\partial_\omega R_2 = -|u|_-\tau_c\f$;
 * - \f$R_3 = n_0 k_0|u|_+ - (d + |u|_+)k_c - \dot k_c\f$:
 *   \f$\partial_{e_c} R_3 = n_0' k_0|u|_+\f$,
 *   \f$\partial_{\dot e_c} R_3 = -\alpha\operatorname{sign}(\dot e_c)k_c\f$,
 *   \f$\partial_{k_c} R_3 = -d - |u|_+\f$,
 *   \f$\partial_{\dot k_c} R_3 = -1\f$,
 *   \f$\partial_\omega R_3 = -|u|_- k_c\f$;
 * - \f$R_4 = m_0 - \omega - \alpha_r\dot\omega\f$:
 *   \f$\partial_{e_c} R_4 = m_0'\f$, \f$\partial_\omega R_4 = -1\f$,
 *   \f$\partial_{\dot\omega} R_4 = -\alpha_r\f$;
 * - \f$R_5 = \mu\dot e_c + \tau_c - \sigma_c\f$:
 *   \f$\partial_{\dot e_c} R_5 = \mu\f$, \f$\partial_{\tau_c} R_5 = 1\f$,
 *   \f$\partial_{\sigma_c} R_5 = -1\f$;
 * - \f$R_6 = \tau g - \sigma_s\f$: \f$\partial_\tau R_6 = g\f$,
 *   \f$\partial_{e_c} R_6 = 2\tau\f$, \f$\partial_{\sigma_s} R_6 = -1\f$;
 * - \f$R_7 = \sigma_s/E_s - e_s\f$: \f$\partial_{\sigma_s} R_7 = 1/E_s\f$,
 *   \f$\partial_{e_s} R_7 = -1\f$.
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
  void linearize(const double* a, const double* adot, double lambda,
                 Linearization& lin) const override;

 private:
  double u_plus_ = 0.0;   ///< positive part of the reaction-rate signal
  double u_minus_ = 0.0;  ///< negative part of the reaction-rate signal
};

/**
 * @brief Bestel-Clement-Sorine active stress of \cite genet23
 *
 * The law of \cite genet23 (Eqs. 29-33, after Chapelle et al. 2012) with the
 * fiber strain \f$e_\text{fib} = \lambda - 1\f$ and a linear series spring.
 * Its contractile element is that of `strain_dependent` (\cite caruel13),
 * whose series element is nonlinear in the Green-Lagrange strain and whose
 * relaxation depends on the length.
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
 * Linearization, with the decay \f$d = |\nu| + \alpha|\dot e_c|\f$:
 * - \f$R_0 = \tau - (\tau_c + \mu\dot e_c)/\lambda\f$:
 *   \f$\partial_\tau R_0 = 1\f$, \f$\partial_{\tau_c} R_0 = -1/\lambda\f$,
 *   \f$\partial_{\dot e_c} R_0 = -\mu/\lambda\f$,
 *   \f$\partial_\lambda R_0 = (\tau_c + \mu\dot e_c)/\lambda^2\f$;
 * - \f$R_1 = \dot e_c - (k_s(\lambda - 1 - e_c) - \tau_c)/\mu\f$:
 *   \f$\partial_{e_c} R_1 = k_s/\mu\f$, \f$\partial_{\tau_c} R_1 = 1/\mu\f$,
 *   \f$\partial_{\dot e_c} R_1 = 1\f$,
 *   \f$\partial_\lambda R_1 = -k_s/\mu\f$;
 * - \f$R_2 = \dot\tau_c + d\,\tau_c - \dot e_c k_c - n_0\sigma_0|\nu|_+\f$:
 *   \f$\partial_{e_c} R_2 = -n_0'\sigma_0|\nu|_+\f$,
 *   \f$\partial_{\dot e_c} R_2 = \alpha\operatorname{sign}(\dot e_c)\tau_c -
 * k_c\f$,
 *   \f$\partial_{\tau_c} R_2 = d\f$, \f$\partial_{\dot\tau_c} R_2 = 1\f$,
 *   \f$\partial_{k_c} R_2 = -\dot e_c\f$;
 * - \f$R_3 = \dot k_c + d\,k_c - n_0 k_0|\nu|_+\f$:
 *   \f$\partial_{e_c} R_3 = -n_0' k_0|\nu|_+\f$,
 *   \f$\partial_{\dot e_c} R_3 = \alpha\operatorname{sign}(\dot e_c)k_c\f$,
 *   \f$\partial_{k_c} R_3 = d\f$, \f$\partial_{\dot k_c} R_3 = 1\f$.
 *
 * Parameters: `sigma_0`, `k_s` (default 1e8), `k_0` (default 260e3), `mu`
 * (default 70), `alpha` (default 0), `n0_e_c`, `n0_values` (optional)
 */
class LinearSeriesSpringActiveStress : public ActiveStress {
 public:
  LinearSeriesSpringActiveStress()
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
  void linearize(const double* a, const double* adot, double lambda,
                 Linearization& lin) const override;

 private:
  double nu_abs_ = 0.0;   ///< |nu(t)|
  double nu_plus_ = 0.0;  ///< |nu(t)|_+
};

#endif  // SVZERODSOLVER_MODEL_ACTIVESTRESS_HPP_
