// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

/**
 * @file ActiveStress.h
 * @brief Active stress models for chamber blocks
 */

#ifndef SVZERODSOLVER_MODEL_ACTIVESTRESS_HPP_
#define SVZERODSOLVER_MODEL_ACTIVESTRESS_HPP_

#include <Eigen/Core>
#include <list>
#include <memory>
#include <string>
#include <vector>

#include "ActivationFunction.h"
#include "Parameter.h"
#include "ParameterHolder.h"
#include "SparseSystem.h"

/**
 * @brief Base class for chamber active stress models
 *
 * An ActiveStress owns the dynamics of a chamber's active stress \f$\tau\f$.
 * Subclasses may introduce additional internal variables and equations
 * beyond ChamberSphere's fixed core state, and write directly into the same
 * global SparseSystem that ChamberSphere writes into, using indices that are
 * a fixed offset into ChamberSphere's own \ref Block::global_var_ids and
 * \ref Block::global_eqn_ids arrays (which are passed down unchanged).
 *
 * ChamberSphere::setup_dofs() always registers, in this order:
 * internal variables `radius, velo, stress, active_stress, volume,
 * passive_stress,
 * <extra_var_names()...>` and equations `momentum, spherical_stress,
 * volume_change, acceleration, mass_conservation, pressure_equality,
 * passive_stress, <this model's equations...>`, so the offsets below are
 * invariant across ActiveStress subclasses.
 */
class ActiveStress : public ParameterHolder {
 public:
  /**
   * @name Fixed layout of ChamberSphere's always-present ("core") state
   *
   * Indices into the full \ref Block::global_var_ids / \ref
   * Block::global_eqn_ids arrays that ChamberSphere passes down unchanged.
   */
  ///@{
  static constexpr int RADIUS_VAR =
      4;  ///< Radius position within ChamberSphere's own global_var_ids array
  static constexpr int VELO_VAR =
      5;  ///< Radial velocity position within ChamberSphere's own
          ///< global_var_ids array
          ///<
  static constexpr int STRESS_VAR =
      6;  ///< Total wall stress position within ChamberSphere's own
          ///< global_var_ids array
          ///<
  static constexpr int TAU_VAR =
      7;  ///< Active stress \f$\tau\f$ position within ChamberSphere's own
          ///< global_var_ids array
          ///<
  static constexpr int VOLUME_VAR =
      8;  ///< Volume position within ChamberSphere's own global_var_ids array
  static constexpr int PASSIVE_STRESS_VAR =
      9;  ///< Passive wall stress position within ChamberSphere's own
          ///< global_var_ids array
          ///<
  static constexpr int NUM_CORE_VARS =
      10;  ///< A subclass's own extra
           ///< variables start at
           ///< global_var_ids[NUM_CORE_VARS + i]
  static constexpr int NUM_CORE_EQNS =
      7;  ///< A subclass's own extra
          ///< equations start at
          ///< global_eqn_ids[NUM_CORE_EQNS + i]
  ///@}

  /**
   * @brief Construct an ActiveStress model
   *
   * @param props Parameter name/spec pairs for this active stress type
   */
  explicit ActiveStress(
      const std::vector<std::pair<std::string, InputParameter>>& props);

  /**
   * @brief Number of internal variables this model owns beyond
   * ChamberSphere's fixed 6 (radius, velo, stress, active_stress,
   * volume, passive_stress)
   *
   * @return Number of extra internal variables
   */
  virtual int num_extra_vars() const { return 0; }

  /**
   * @brief Names of this model's extra internal variables, in registration
   * order (appended immediately after "passive_stress")
   *
   * @return List of extra internal variable names
   */
  virtual std::list<std::string> extra_var_names() const { return {}; }

  /**
   * @brief Number of equations this model owns beyond ChamberSphere's fixed
   * 7 core equations. Every ActiveStress supplies at least the equation
   * defining \f$\tau\f$.
   *
   * @return Number of extra equations
   */
  virtual int num_extra_equations() const = 0;

  /**
   * @brief Constant (time- and solution-independent) contributions
   *
   * @param system System to update contributions at
   * @param global_eqn_ids ChamberSphere's global equation indices
   * @param global_var_ids ChamberSphere's global variable indices
   */
  virtual void update_constant(SparseSystem& system,
                               const std::vector<int>& global_eqn_ids,
                               const std::vector<int>& global_var_ids) {
    (void)system;
    (void)global_eqn_ids;
    (void)global_var_ids;
  }

  /**
   * @brief Time-dependent contributions
   *
   * @param system System to update contributions at
   * @param global_eqn_ids ChamberSphere's global equation indices
   * @param global_var_ids ChamberSphere's global variable indices
   * @param activation_signal The chamber's ActivationFunction evaluated at
   * the current time
   */
  virtual void update_time(SparseSystem& system,
                           const std::vector<int>& global_eqn_ids,
                           const std::vector<int>& global_var_ids,
                           double activation_signal) {
    (void)system;
    (void)global_eqn_ids;
    (void)global_var_ids;
    (void)activation_signal;
  }

  /**
   * @brief Solution-dependent contributions
   *
   * @param system System to update contributions at
   * @param global_eqn_ids ChamberSphere's global equation indices
   * @param global_var_ids ChamberSphere's global variable indices
   * @param y Current solution
   * @param dy Current derivative of the solution
   * @param radius0 Reference radius \f$r_0\f$ of the chamber
   * @param activation_signal The chamber's ActivationFunction evaluated at
   * the current time
   */
  virtual void update_solution(
      SparseSystem& system, const std::vector<int>& global_eqn_ids,
      const std::vector<int>& global_var_ids,
      const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
      const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy, double radius0,
      double activation_signal) = 0;

  /**
   * @brief Type of activation signal this active stress model needs
   *
   * Checked against the chamber's ActivationFunction::signal_type in
   * ChamberSphere::setup_dofs to catch mismatched pairings (e.g.
   * strain_dependent active stress with a normalized activation function).
   *
   * @return The required ActivationSignalType
   */
  virtual ActivationSignalType required_signal_type() const = 0;

  /**
   * @brief Factory: create an active stress model from a type string
   *
   * @param type_str One of: "strain_independent", "strain_dependent"
   * @return Unique pointer to the created active stress model
   */
  static std::unique_ptr<ActiveStress> create(const std::string& type_str);
};

/**
 * @brief Strain-independent active stress model
 *
 * Owns 1 extra equation (the \f$\tau\f$ ODE) and 0 extra variables:
 * \f[
 * \dot{\tau} + |a| \tau - \sigma_\text{max} a_+ = 0, \quad a_+ = \max(a, 0),
 * \quad a = f \alpha_\text{max} + (1 - f) \alpha_\text{min}
 * \f]
 * where \f$f \in [0, 1]\f$ is the activation function.
 *
 * Parameters: `alpha_max`, `alpha_min`, `sigma_max`
 */
class StrainIndependentActiveStress : public ActiveStress {
 public:
  StrainIndependentActiveStress()
      : ActiveStress({{"alpha_max", InputParameter()},
                      {"alpha_min", InputParameter()},
                      {"sigma_max", InputParameter()}}) {}

  int num_extra_equations() const override { return 1; }

  ActivationSignalType required_signal_type() const override {
    return ActivationSignalType::normalized;
  }

  void update_constant(SparseSystem& system,
                       const std::vector<int>& global_eqn_ids,
                       const std::vector<int>& global_var_ids) override;

  void update_time(SparseSystem& system, const std::vector<int>& global_eqn_ids,
                   const std::vector<int>& global_var_ids,
                   double activation_signal) override;

  void update_solution(SparseSystem& system,
                       const std::vector<int>& global_eqn_ids,
                       const std::vector<int>& global_var_ids,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy,
                       double radius0, double activation_signal) override;

 private:
  double act_ = 0.0;       ///< Rate coefficient
  double act_plus_ = 0.0;  ///< Positive part of act_, act_plus_ = max(act_, 0)
};

/**
 * @brief Strain-dependent active stress model from \cite caruel13
 *
 * Owns 8 extra equations and 7 extra variables (\f$e_c, \tau_c, k_c,
 * \omega, \sigma_c, \sigma_s, e_s\f$, in that order):
 * \f[
 * \tau = E_s \frac{e_{1D} - e_c}{(1 + 2 e_c)^2}
 * \f]
 * with \f$e_{1D}\f$ the strain in fiber direction (a function of the
 * chamber's radius), and evolution equations
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

  int num_extra_vars() const override { return 7; }
  std::list<std::string> extra_var_names() const override {
    return {"e_c", "tau_c", "k_c", "omega", "sigma_c", "sigma_s", "e_s"};
  }
  int num_extra_equations() const override { return 8; }

  ActivationSignalType required_signal_type() const override {
    return ActivationSignalType::rate;
  }

  void update_constant(SparseSystem& system,
                       const std::vector<int>& global_eqn_ids,
                       const std::vector<int>& global_var_ids) override;

  void update_time(SparseSystem& system, const std::vector<int>& global_eqn_ids,
                   const std::vector<int>& global_var_ids,
                   double activation_signal) override;

  void update_solution(SparseSystem& system,
                       const std::vector<int>& global_eqn_ids,
                       const std::vector<int>& global_var_ids,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy,
                       double radius0, double activation_signal) override;

 private:
  /**
   * @name Extra variable/equation positions within the full ChamberSphere
   * global_var_ids/global_eqn_ids arrays
   *
   * Reproduce the exact global_var_ids/global_eqn_ids indices that the
   * ChamberSphere used.
   */
  ///@{
  static constexpr int E_C_VAR = NUM_CORE_VARS + 0;      // 10
  static constexpr int TAU_C_VAR = NUM_CORE_VARS + 1;    // 11
  static constexpr int K_C_VAR = NUM_CORE_VARS + 2;      // 12
  static constexpr int OMEGA_VAR = NUM_CORE_VARS + 3;    // 13
  static constexpr int SIGMA_C_VAR = NUM_CORE_VARS + 4;  // 14
  static constexpr int SIGMA_S_VAR = NUM_CORE_VARS + 5;  // 15
  static constexpr int E_S_VAR = NUM_CORE_VARS + 6;      // 16
  static constexpr int TAU_EQN = NUM_CORE_EQNS + 0;      // 7
  static constexpr int OMEGA_EQN = NUM_CORE_EQNS + 1;    // 8
  static constexpr int E_C_EQN = NUM_CORE_EQNS + 2;      // 9
  static constexpr int K_C_EQN = NUM_CORE_EQNS + 3;      // 10
  static constexpr int TAU_C_EQN = NUM_CORE_EQNS + 4;    // 11
  static constexpr int SIGMA_C_EQN = NUM_CORE_EQNS + 5;  // 12
  static constexpr int SIGMA_S_EQN = NUM_CORE_EQNS + 6;  // 13
  static constexpr int E_S_EQN = NUM_CORE_EQNS + 7;      // 14
  ///@}

  /**
   * @brief Update the strain-dependent coefficients used across this
   * model's equations
   *
   * @param e_c Current sarcomere element strain
   */
  void update_active_stress_values(double e_c);

  double n_0_ = 0.0;        // activation strain-dependence
  double m_0_ = 0.0;        // relaxation strain-dependence
  double dn_0_de_c_ = 0.0;  // d(n_0)/d(e_c), local slope of n0_e_c/n0_values
  double dm_0_de_c_ = 0.0;  // d(m_0)/d(e_c), local slope of m0_e_c/m0_values
  double u_plus_ = 0.0;     // positive part of the reaction-rate signal
  double u_minus_ = 0.0;    // negative part of the reaction-rate signal
};

#endif  // SVZERODSOLVER_MODEL_ACTIVESTRESS_HPP_
