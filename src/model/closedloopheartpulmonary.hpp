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
/**
 * @file closedloopheartpulmonary.hpp
 * @brief MODEL::ClosedLoopHeartPulmonary source file
 */
#ifndef SVZERODSOLVER_MODEL_CLOSEDLOOPHEARTPULMONARY_HPP_
#define SVZERODSOLVER_MODEL_CLOSEDLOOPHEARTPULMONARY_HPP_

#include <math.h>

#include "../algebra/sparsesystem.hpp"
#include "../algebra/state.hpp"
#include "block.hpp"
#define PI 3.14159265

namespace MODEL {
/**
 * @brief Heart and pulmonary circulation model
 *
 * Models the mechanics of the 4 heart chambers and pulmonary circulation
 *
 * Reference for equations and model structure: Sankaran, S., Moghadam, M. E.,
 * Kahn, A. M., Tseng, E. E., Guccione, J. M., & Marsden, A. L. (2012).
 * Patient-specific multiscale modeling of blood flow for coronary artery bypass
 * graft surgery. Annals of Biomedical Engineering, 40(10), 2228–2242.
 * https://doi.org/10.1007/s10439-012-0579-3
 *
 * TODO: Equations and circuit diagram
 *
 * @tparam T Scalar type (e.g. `float`, `double`)
 */
template <typename T>
class ClosedLoopHeartPulmonary : public Block<T> {
 public:
  /**
   * @brief Parameters of the element.
   *
   * Struct containing all constant and/or time-dependent parameters of the
   * element.
   */
  struct Parameters : public Block<T>::Parameters {
    T Tsa;      // Fractions of cardiac cycle (not sure)
    T tpwave;   // Fraction of cardiac cycle (P-wave)
    T Erv_s;    // Scaling for right ventricle elastance
    T Elv_s;    // Scaling for left ventricle elastance
    T iml;      // Scaling for intramyocardial pressure (left coronaries)
    T imr;      // Scaling for intramyocardial pressure (right coronaries)
    T Lra_v;    // Right atrium inductance
    T Rra_v;    // Right atrium outflow resistance
    T Lrv_a;    // Right ventricle inductance
    T Rrv_a;    // Right ventricle outflow resistance
    T Lla_v;    // Left atrium inductance
    T Rla_v;    // Left atrium outflow resistance
    T Llv_a;    // Left ventricle inductance
    T Rlv_ao;   // Left ventricle outflow resistance
    T Vrv_u;    // Right ventricle unstressed volume
    T Vlv_u;    // Right ventricle unstressed volume
    T Rpd;      // Pulmonary resistance
    T Cp;       // Pulmonary capacitance
    T Cpa;      // Aortic capacitance
    T Kxp_ra;   // Right atrium pressure-volume relationship (?)
    T Kxv_ra;   // Right atrium pressure-volume relationship (?)
    T Kxp_la;   // Left atrium pressure-volume relationship (?)
    T Kxv_la;   // Left atrium pressure-volume relationship (?)
    T Emax_ra;  // Right atrium elastance (?)
    T Emax_la;  // Left atrium elastance (?)
    T Vaso_ra;  // Right atrium rest volume (?)
    T Vaso_la;  // Left atrium rest volume (?)
  };

  /**
   * @brief Construct a new ClosedLoopHeartPulmonary object
   *
   * @param heart_parameters Map/dictionary containing the 27 heart parameters
   * @param name Name
   */
  // ClosedLoopHeartPulmonary(std::map<std::string, T> heart_parameters,
  // std::string name);
  ClosedLoopHeartPulmonary(std::map<std::string, T> heart_parameters,
                           T cycle_period, std::string name);

  /**
   * @brief Destroy the ClosedLoopHeartPulmonary object
   *
   */
  ~ClosedLoopHeartPulmonary();

  /**
   * @brief Set up the degrees of freedom (DOF) of the block
   *
   * Set \ref global_var_ids and \ref global_eqn_ids of the element based on the
   * number of equations and the number of internal variables of the
   * element.
   *
   * @param dofhandler Degree-of-freedom handler to register variables and
   * equations at
   */
  void setup_dofs(DOFHandler &dofhandler);

  /**
   * @brief Update for model-dependent variables that are set at the end of
   * model construction
   *
   * @param model Model object to access model-dependent variables
   */
  // void update_model_dependent_params(MODEL::Model<T> &model);

  /**
   * @brief Return parameter values
   *
   * @param message String to identify different requests
   */
  void get_parameter_value(std::string message, T &param_value);

  /**
   * @brief Set block-specific initial conditions
   *
   * @param state State vector containing y and ydot
   */
  void set_ICs(ALGEBRA::State<T> &state);

  /**
   * @brief Update the constant contributions of the element in a sparse system
   *
   * @param system System to update contributions at
   */
  void update_constant(ALGEBRA::SparseSystem<T> &system);

  /**
   * @brief Update the time-dependent contributions of the element in a sparse
   * system
   *
   * @param system System to update contributions at
   * @param time Current time
   */
  void update_time(ALGEBRA::SparseSystem<T> &system, T time);

  /**
   * @brief Update the solution-dependent contributions of the element in a
   * sparse system
   *
   * @param system System to update contributions at
   * @param y Current solution
   */
  void update_solution(ALGEBRA::SparseSystem<T> &system,
                       Eigen::Matrix<T, Eigen::Dynamic, 1> &y);

  /**
   * @brief Number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  std::map<std::string, int> num_triplets = {
      {"F", 33},
      {"E", 10},
      {"D", 2},
  };

  /**
   * @brief Get number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  std::map<std::string, int> get_num_triplets();

  /**
   * @brief Update the atrial activation and LV/RV elastance functions which
   * depend on time
   *
   * @param time Current time
   */
  void get_activation_and_elastance_functions(T time);

  /**
   * @brief Compute sub-expressions that are part of atrial elastance and
   * depends on atrial volume from the solution vector
   *
   * @param y Current solution
   */
  void get_psi_ra_la(Eigen::Matrix<T, Eigen::Dynamic, 1> &y);

  /**
   * @brief Valve positions for each heart chamber
   *
   * @param y Current solution
   * @param valves Vector containing valve positions
   */
  void get_valve_positions(Eigen::Matrix<T, Eigen::Dynamic, 1> &y);

  T cardiac_cycle_period;

 private:
  Parameters params;
  // Below variables change every timestep and are then combined with
  // expressions that are updated with solution
  T AA;   // Atrial activation function
  T Elv;  // LV elastance
  T Erv;  // RV elastance
  T psi_ra, psi_la, psi_ra_derivative,
      psi_la_derivative;  // Expressions for atrial activation
  T valves[16];
};

template <typename T>
ClosedLoopHeartPulmonary<T>::ClosedLoopHeartPulmonary(
    std::map<std::string, T> heart_parameters, T cycle_period, std::string name)
    : Block<T>(name) {
  this->name = name;
  this->cardiac_cycle_period = cycle_period;
  this->params.Tsa = heart_parameters["Tsa"];
  this->params.tpwave = heart_parameters["tpwave"];
  this->params.Erv_s = heart_parameters["Erv_s"];
  this->params.Elv_s = heart_parameters["Elv_s"];
  this->params.iml = heart_parameters["iml"];
  this->params.imr = heart_parameters["imr"];
  this->params.Lra_v = heart_parameters["Lra_v"];
  this->params.Rra_v = heart_parameters["Rra_v"];
  this->params.Lrv_a = heart_parameters["Lrv_a"];
  this->params.Rrv_a = heart_parameters["Rrv_a"];
  this->params.Lla_v = heart_parameters["Lla_v"];
  this->params.Rla_v = heart_parameters["Rla_v"];
  this->params.Llv_a = heart_parameters["Llv_a"];
  this->params.Rlv_ao = heart_parameters["Rlv_ao"];
  this->params.Vrv_u = heart_parameters["Vrv_u"];
  this->params.Vlv_u = heart_parameters["Vlv_u"];
  this->params.Rpd = heart_parameters["Rpd"];
  this->params.Cp = heart_parameters["Cp"];
  this->params.Cpa = heart_parameters["Cpa"];
  this->params.Kxp_ra = heart_parameters["Kxp_ra"];
  this->params.Kxv_ra = heart_parameters["Kxv_ra"];
  this->params.Kxp_la = heart_parameters["Kxp_la"];
  this->params.Kxv_la = heart_parameters["Kxv_la"];
  this->params.Emax_ra = heart_parameters["Emax_ra"];
  this->params.Emax_la = heart_parameters["Emax_la"];
  this->params.Vaso_ra = heart_parameters["Vaso_ra"];
  this->params.Vaso_la = heart_parameters["Vaso_la"];
}

template <typename T>
ClosedLoopHeartPulmonary<T>::~ClosedLoopHeartPulmonary() {}

template <typename T>
void ClosedLoopHeartPulmonary<T>::setup_dofs(DOFHandler &dofhandler) {
  // Block<T>::setup_dofs_(dofhandler, 14, 12);
  Block<T>::setup_dofs_(dofhandler, 14,
                        {"V_RA", "Q_RA", "P_RV", "V_RV", "Q_RV", "P_pul",
                         "P_LA", "V_LA", "Q_LA", "P_LV", "V_LV", "Q_LV"});
}

template <typename T>
void ClosedLoopHeartPulmonary<T>::update_constant(
    ALGEBRA::SparseSystem<T> &system) {
  // DOF 2, Eq 1: Aortic pressure
  system.E.coeffRef(this->global_eqn_ids[1], this->global_var_ids[2]) =
      params.Cpa;
  // DOF 4, Eq 2: Right atrium volume
  system.E.coeffRef(this->global_eqn_ids[2], this->global_var_ids[4]) = 1.0;
  // DOF 5, Eq 3: Right atrium outflow
  system.E.coeffRef(this->global_eqn_ids[3], this->global_var_ids[5]) =
      params.Lra_v;
  // DOF 7, Eq 5: Right ventricle volume
  system.E.coeffRef(this->global_eqn_ids[5], this->global_var_ids[7]) = 1.0;
  // DOF 8, Eq 6: Right ventricle outflow
  system.E.coeffRef(this->global_eqn_ids[6], this->global_var_ids[8]) =
      params.Lrv_a;
  // DOF 9, Eq 7: Pulmonary pressure
  system.E.coeffRef(this->global_eqn_ids[7], this->global_var_ids[9]) =
      params.Cp;
  // DOF 11, Eq 9: Left atrium volume
  system.E.coeffRef(this->global_eqn_ids[9], this->global_var_ids[11]) = 1.0;
  // DOF 12, Eq 10: Left atrium outflow
  system.E.coeffRef(this->global_eqn_ids[10], this->global_var_ids[12]) =
      params.Lla_v;
  // DOF 14, Eq 12: Left ventricle volume
  system.E.coeffRef(this->global_eqn_ids[12], this->global_var_ids[14]) = 1.0;
  // DOF 15, Eq 13: Left ventricle outflow
  system.E.coeffRef(this->global_eqn_ids[13], this->global_var_ids[15]) =
      params.Llv_a;
}

template <typename T>
void ClosedLoopHeartPulmonary<T>::update_time(ALGEBRA::SparseSystem<T> &system,
                                              T time) {
  this->get_activation_and_elastance_functions(time);
}

template <typename T>
void ClosedLoopHeartPulmonary<T>::update_solution(
    ALGEBRA::SparseSystem<T> &system, Eigen::Matrix<T, Eigen::Dynamic, 1> &y) {
  this->get_psi_ra_la(y);
  this->get_valve_positions(y);

  // F and C matrices depend on time and solution
  // Specifying all terms here, including constant terms (which can instead be
  // specified in update_constant) for readability (Doesn't seem to make a
  // difference to compute time) DOF IDs are arranged as inflow
  // [P_in,Q_in,P_out,Q_out,internal variables...]

  // DOF 0, Eq 0: Right atrium pressure
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[0]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[4]) =
      -this->AA * params.Emax_ra;
  system.C(this->global_eqn_ids[0]) =
      this->AA * params.Emax_ra * params.Vaso_ra + psi_ra * (this->AA - 1.0);
  system.D.coeffRef(this->global_eqn_ids[0], this->global_var_ids[4]) =
      psi_ra_derivative * (this->AA - 1.0);

  // DOF 1: Flow into right atrium (no equation)

  // DOF 2, Eq 1: Aortic pressure
  system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[15]) =
      -valves[15];
  system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[3]) = 1.0;

  // DOF 3: Flow into aorta (no equation)

  // DOF 4, Eq 2: Right atrium volume
  system.F.coeffRef(this->global_eqn_ids[2], this->global_var_ids[5]) =
      1.0 * valves[5];
  system.F.coeffRef(this->global_eqn_ids[2], this->global_var_ids[1]) = -1.0;

  // DOF 5, Eq 3: Right atrium outflow
  system.F.coeffRef(this->global_eqn_ids[3], this->global_var_ids[5]) =
      params.Rra_v * valves[5];
  system.F.coeffRef(this->global_eqn_ids[3], this->global_var_ids[0]) = -1.0;
  system.F.coeffRef(this->global_eqn_ids[3], this->global_var_ids[6]) = 1.0;

  // DOF 6, Eq 4: Right ventricle pressure
  system.F.coeffRef(this->global_eqn_ids[4], this->global_var_ids[6]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[4], this->global_var_ids[7]) =
      -this->Erv;
  system.C(this->global_eqn_ids[4]) = this->Erv * params.Vrv_u;

  // DOF 7, Eq 5: Right ventricle volume
  system.F.coeffRef(this->global_eqn_ids[5], this->global_var_ids[5]) =
      -1.0 * valves[5];
  system.F.coeffRef(this->global_eqn_ids[5], this->global_var_ids[8]) =
      1.0 * valves[8];

  // DOF 8, Eq 6: Right ventricle outflow
  system.F.coeffRef(this->global_eqn_ids[6], this->global_var_ids[6]) = -1.0;
  system.F.coeffRef(this->global_eqn_ids[6], this->global_var_ids[9]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[6], this->global_var_ids[8]) =
      params.Rrv_a * valves[8];

  // DOF 9, Eq 7: Pulmonary pressure
  system.F.coeffRef(this->global_eqn_ids[7], this->global_var_ids[8]) =
      -valves[8];
  system.F.coeffRef(this->global_eqn_ids[7], this->global_var_ids[9]) =
      1.0 / params.Rpd;
  system.F.coeffRef(this->global_eqn_ids[7], this->global_var_ids[10]) =
      -1.0 / params.Rpd;

  // DOF 10, Eq 8: Left atrium pressure
  system.F.coeffRef(this->global_eqn_ids[8], this->global_var_ids[10]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[8], this->global_var_ids[11]) =
      -this->AA * params.Emax_la;
  system.C(this->global_eqn_ids[8]) =
      this->AA * params.Emax_la * params.Vaso_la + psi_la * (this->AA - 1.0);
  system.D.coeffRef(this->global_eqn_ids[8], this->global_var_ids[11]) =
      psi_la_derivative * (this->AA - 1.0);

  // DOF 11, Eq 9: Left atrium volume
  system.F.coeffRef(this->global_eqn_ids[9], this->global_var_ids[8]) =
      -1.0 * valves[8];
  system.F.coeffRef(this->global_eqn_ids[9], this->global_var_ids[12]) =
      1.0 * valves[12];

  // DOF 12, Eq 10: Left atrium outflow
  system.F.coeffRef(this->global_eqn_ids[10], this->global_var_ids[10]) = -1.0;
  system.F.coeffRef(this->global_eqn_ids[10], this->global_var_ids[13]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[10], this->global_var_ids[12]) =
      params.Rla_v * valves[12];

  // DOF 13, Eq 11: Left ventricle pressure
  system.F.coeffRef(this->global_eqn_ids[11], this->global_var_ids[13]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[11], this->global_var_ids[14]) =
      -this->Elv;
  system.C(this->global_eqn_ids[11]) = this->Elv * params.Vlv_u;

  // DOF 14, Eq 12: Left ventricle volume
  system.F.coeffRef(this->global_eqn_ids[12], this->global_var_ids[12]) =
      -1.0 * valves[12];
  system.F.coeffRef(this->global_eqn_ids[12], this->global_var_ids[15]) =
      1.0 * valves[15];

  // DOF 15, Eq 13: Left ventricle outflow
  system.F.coeffRef(this->global_eqn_ids[13], this->global_var_ids[13]) = -1.0;
  system.F.coeffRef(this->global_eqn_ids[13], this->global_var_ids[2]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[13], this->global_var_ids[15]) =
      params.Rlv_ao * valves[15];
}

template <typename T>
void ClosedLoopHeartPulmonary<T>::get_activation_and_elastance_functions(
    T time) {
  T T_cardiac = this->cardiac_cycle_period;
  T Tsa = T_cardiac * params.Tsa;
  T tpwave = T_cardiac / params.tpwave;
  T t_in_cycle = fmod(time, T_cardiac);

  // Activation function
  AA = 0.0;
  if (t_in_cycle <= tpwave) {
    AA = (0.5) * (1.0 - cos(2.0 * PI * (t_in_cycle - tpwave + Tsa) / Tsa));
  } else if ((t_in_cycle >= (T_cardiac - Tsa) + tpwave) and
             (t_in_cycle < T_cardiac)) {
    AA =
        (0.5) *
        (1.0 - cos(2.0 * PI * (t_in_cycle - tpwave - (T_cardiac - Tsa)) / Tsa));
  } else {
    AA = 0.0;
  }

  // Elastance modes (copied from J. Tran's tuning framework)
  const int num_elast_modes = 25;
  T Ft_elastance[num_elast_modes][2] = {
      {0.283748803, 0.000000000},   {0.031830626, -0.374299825},
      {-0.209472400, -0.018127770}, {0.020520047, 0.073971113},
      {0.008316883, -0.047249597},  {-0.041677660, 0.003212163},
      {0.000867323, 0.019441411},   {-0.001675379, -0.005565534},
      {-0.011252277, 0.003401432},  {-0.000414677, 0.008376795},
      {0.000253749, -0.000071880},  {-0.002584966, 0.001566861},
      {0.000584752, 0.003143555},   {0.000028502, -0.000024787},
      {0.000022961, -0.000007476},  {0.000018735, -0.000001281},
      {0.000015573, 0.000001781},   {0.000013133, 0.000003494},
      {0.000011199, 0.000004507},   {0.000009634, 0.000005117},
      {0.000008343, 0.000005481},   {0.000007265, 0.000005687},
      {0.000006354, 0.000005789},   {0.000005575, 0.000005821},
      {0.000004903, 0.000005805}};

  // RV and LV elastance
  T Elv_i = 0.0;
  for (auto i = 0; i < num_elast_modes; i++)
    Elv_i = Elv_i +
            (Ft_elastance[i][0]) * cos(2.0 * PI * i * t_in_cycle / T_cardiac) -
            (Ft_elastance[i][1]) * sin(2.0 * PI * i * t_in_cycle / T_cardiac);

  Elv = Elv_i * params.Elv_s;
  Erv = Elv_i * params.Erv_s;
}

template <typename T>
void ClosedLoopHeartPulmonary<T>::get_psi_ra_la(
    Eigen::Matrix<T, Eigen::Dynamic, 1> &y) {
  auto RA_volume = y[this->global_var_ids[4]];
  auto LA_volume = y[this->global_var_ids[11]];
  psi_ra =
      params.Kxp_ra * (exp((RA_volume - params.Vaso_ra) * params.Kxv_ra) - 1.0);
  psi_la =
      params.Kxp_la * (exp((LA_volume - params.Vaso_la) * params.Kxv_la) - 1.0);

  psi_ra_derivative = params.Kxp_ra *
                      exp((RA_volume - params.Vaso_ra) * params.Kxv_ra) *
                      params.Kxv_ra;
  psi_la_derivative = params.Kxp_la *
                      exp((LA_volume - params.Vaso_la) * params.Kxv_la) *
                      params.Kxv_la;
}

template <typename T>
void ClosedLoopHeartPulmonary<T>::get_valve_positions(
    Eigen::Matrix<T, Eigen::Dynamic, 1> &y) {
  std::fill(valves, valves + 16, 1.0);

  // RA to RV
  auto pressure_ra = y[this->global_var_ids[0]];
  auto pressure_rv = y[this->global_var_ids[6]];
  auto outflow_ra = y[this->global_var_ids[5]];
  if ((pressure_ra <= pressure_rv) and (outflow_ra <= 0.0)) {
    valves[5] = 0.0;
    y[this->global_var_ids[5]] = 0.0;
  }

  // RV to pulmonary
  auto pressure_pulmonary = y[this->global_var_ids[9]];
  auto outflow_rv = y[this->global_var_ids[8]];
  if ((pressure_rv <= pressure_pulmonary) and (outflow_rv <= 0.0)) {
    valves[8] = 0.0;
    y[this->global_var_ids[8]] = 0.0;
  }

  // LA to LV
  auto pressure_la = y[this->global_var_ids[10]];
  auto pressure_lv = y[this->global_var_ids[13]];
  auto outflow_la = y[this->global_var_ids[12]];
  if ((pressure_la <= pressure_lv) and (outflow_la <= 0.0)) {
    valves[12] = 0.0;
    y[this->global_var_ids[12]] = 0.0;
  }

  // LV to aorta
  auto pressure_aorta = y[this->global_var_ids[2]];
  auto outflow_lv = y[this->global_var_ids[15]];
  if ((pressure_lv <= pressure_aorta) and (outflow_lv <= 0.0)) {
    valves[15] = 0.0;
    y[this->global_var_ids[15]] = 0.0;
  }
}

template <typename T>
void ClosedLoopHeartPulmonary<T>::set_ICs(ALGEBRA::State<T> &state) {
  state.y[this->global_var_ids[4]] = 38.43;   // RA vol.
  state.y[this->global_var_ids[7]] = 96.07;   // RV vol.
  state.y[this->global_var_ids[11]] = 38.43;  // LA vol.
  state.y[this->global_var_ids[14]] = 96.07;  // LV vol.
  state.y[this->global_var_ids[9]] = 8.0;     // Pulm pressure
  // Below ICs likely are not needed (but retained as comments in case they are)
  // state.y[this->global_var_ids[0]] = 4.72;   // RA pressure
  // state.y[this->global_var_ids[6]] = 14.58;  // RV pressure
  // state.y[this->global_var_ids[10]] = 6.09;  // LA pressure
  // state.y[this->global_var_ids[13]] = 22.39; // LV pressure
}

template <typename T>
void ClosedLoopHeartPulmonary<T>::get_parameter_value(std::string message,
                                                      T &param_value) {
  if (message == "iml") {
    param_value = params.iml;
  } else if (message == "imr") {
    param_value = params.imr;
  }
}

template <typename T>
std::map<std::string, int> ClosedLoopHeartPulmonary<T>::get_num_triplets() {
  return num_triplets;
}

}  // namespace MODEL

#endif  // SVZERODSOLVER_MODEL_CLOSEDLOOPHEARTPULMONARY_HPP_
