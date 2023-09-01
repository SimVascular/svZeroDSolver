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
 * ### Parameters
 *
 * Parameter sequence for constructing this block
 *
 * * `0` Atrial systole time fraction
 * * `1` Time for P-wave
 * * `2` Scaling for right ventricle elastance
 * * `3` Scaling for left ventricle elastance
 * * `4` Scaling for intramyocardial pressure (left coronaries)
 * * `5` Scaling for intramyocardial pressure (right coronaries)
 * * `6` Right atrium inductance
 * * `7` Right atrium outflow resistance
 * * `8` Right ventricle inductance
 * * `9` Right ventricle outflow resistance
 * * `10` Left atrium inductance
 * * `11` Left atrium outflow resistance
 * * `12` Left ventricle inductance
 * * `13` Left ventricle outflow resistance
 * * `14` Right ventricle unstressed volume
 * * `15` Left ventricle unstressed volume
 * * `16` Pulmonary resistance
 * * `17` Pulmonary capacitance
 * * `18` Aortic capacitance
 * * `19` Right atrium pressure scaling
 * * `20` Right atrium volume scaling
 * * `21` Left atrium pressure scaling
 * * `22` Left atrium volume scaling
 * * `23` Right atrium elastance
 * * `24` Left atrium elastance
 * * `25` Right atrium resting volume
 * * `26` Left atrium resting volume
 *
 * @tparam T Scalar type (e.g. `float`, `double`)
 */
template <typename T>
class ClosedLoopHeartPulmonary : public Block<T> {
 public:
  // Inherit constructors
  using Block<T>::Block;

  /**
   * @brief Local IDs of the parameters
   *
   */
  enum ParamId {
    TSA = 0,       ///< Fractions of cardiac cycle (not sure)
    TPWAVE = 1,    ///< Fraction of cardiac cycle (P-wave)
    ERV_S = 2,     ///< Scaling for right ventricle elastance
    ELV_S = 3,     ///< Scaling for left ventricle elastance
    IML = 4,       ///< Scaling for intramyocardial pressure (left coronaries)
    IMR = 5,       ///< Scaling for intramyocardial pressure (right coronaries)
    LRA_V = 6,     ///< Right atrium inductance
    RRA_V = 7,     ///< Right atrium outflow resistance
    LRV_A = 8,     ///< Right ventricle inductance
    RRV_A = 9,     ///< Right ventricle outflow resistance
    LLA_V = 10,    ///< Left atrium inductance
    RLA_V = 11,    ///< Left atrium outflow resistance
    LLV_A = 12,    ///< Left ventricle inductance
    RLV_AO = 13,   ///< Left ventricle outflow resistance
    VRV_U = 14,    ///< Right ventricle unstressed volume
    VLV_U = 15,    ///< Left ventricle unstressed volume
    RPD = 16,      ///< Pulmonary resistance
    CP = 17,       ///< Pulmonary capacitance
    CPA = 18,      ///< Aortic capacitance
    KXP_RA = 19,   ///< Right atrium pressure-volume relationship (?)
    KXV_RA = 20,   ///< Right atrium pressure-volume relationship (?)
    KXP_LA = 21,   ///< Left atrium pressure-volume relationship (?)
    KXV_LA = 22,   ///< Left atrium pressure-volume relationship (?)
    EMAX_RA = 23,  ///< Right atrium elastance (?)
    EMAX_LA = 24,  ///< Left atrium elastance (?)
    VASO_RA = 25,  ///< Right atrium rest volume (?)
    VASO_LA = 26,  ///< Left atrium rest volume (?)
  };

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
   * @brief Update the constant contributions of the element in a sparse
   system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  void update_constant(ALGEBRA::SparseSystem<T> &system,
                       std::vector<T> &parameters);

  /**
   * @brief Update the time-dependent contributions of the element in a sparse
   * system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  void update_time(ALGEBRA::SparseSystem<T> &system,
                   std::vector<T> &parameters);

  /**
   * @brief Update the solution-dependent contributions of the element in a
   * sparse system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   * @param y Current solution
   * @param dy Current derivate of the solution
   */
  void update_solution(ALGEBRA::SparseSystem<T> &system,
                       std::vector<T> &parameters,
                       Eigen::Matrix<T, Eigen::Dynamic, 1> &y,
                       Eigen::Matrix<T, Eigen::Dynamic, 1> &dy);

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

 private:
  // Below variables change every timestep and are then combined with
  // expressions that are updated with solution
  T AA;   // Atrial activation function
  T Elv;  // LV elastance
  T Erv;  // RV elastance
  T psi_ra, psi_la, psi_ra_derivative,
      psi_la_derivative;  // Expressions for atrial activation
  T valves[16];

  /**
   * @brief Update the atrial activation and LV/RV elastance functions which
   * depend on time
   *
   * @param parameters Parameters of the model
   */
  void get_activation_and_elastance_functions(std::vector<T> &parameters);

  /**
   * @brief Compute sub-expressions that are part of atrial elastance and
   * depends on atrial volume from the solution vector
   *
   * @param parameters Parameters of the model
   * @param y Current solution
   */
  void get_psi_ra_la(std::vector<T> &parameters,
                     Eigen::Matrix<T, Eigen::Dynamic, 1> &y);

  /**
   * @brief Valve positions for each heart chamber
   *
   * @param y Current solution
   */
  void get_valve_positions(Eigen::Matrix<T, Eigen::Dynamic, 1> &y);
};

template <typename T>
void ClosedLoopHeartPulmonary<T>::setup_dofs(DOFHandler &dofhandler) {
  Block<T>::setup_dofs_(dofhandler, 14,
                        {"V_RA", "Q_RA", "P_RV", "V_RV", "Q_RV", "P_pul",
                         "P_LA", "V_LA", "Q_LA", "P_LV", "V_LV", "Q_LV"});
}

template <typename T>
void ClosedLoopHeartPulmonary<T>::update_constant(
    ALGEBRA::SparseSystem<T> &system, std::vector<T> &parameters) {
  // DOF 2, Eq 1: Aortic pressure
  system.E.coeffRef(this->global_eqn_ids[1], this->global_var_ids[2]) =
      parameters[this->global_param_ids[ParamId::CPA]];
  // DOF 4, Eq 2: Right atrium volume
  system.E.coeffRef(this->global_eqn_ids[2], this->global_var_ids[4]) = 1.0;
  // DOF 5, Eq 3: Right atrium outflow
  system.E.coeffRef(this->global_eqn_ids[3], this->global_var_ids[5]) =
      parameters[this->global_param_ids[ParamId::LRA_V]];
  // DOF 7, Eq 5: Right ventricle volume
  system.E.coeffRef(this->global_eqn_ids[5], this->global_var_ids[7]) = 1.0;
  // DOF 8, Eq 6: Right ventricle outflow
  system.E.coeffRef(this->global_eqn_ids[6], this->global_var_ids[8]) =
      parameters[this->global_param_ids[ParamId::LRV_A]];
  // DOF 9, Eq 7: Pulmonary pressure
  system.E.coeffRef(this->global_eqn_ids[7], this->global_var_ids[9]) =
      parameters[this->global_param_ids[ParamId::CP]];
  // DOF 11, Eq 9: Left atrium volume
  system.E.coeffRef(this->global_eqn_ids[9], this->global_var_ids[11]) = 1.0;
  // DOF 12, Eq 10: Left atrium outflow
  system.E.coeffRef(this->global_eqn_ids[10], this->global_var_ids[12]) =
      parameters[this->global_param_ids[ParamId::LLA_V]];
  // DOF 14, Eq 12: Left ventricle volume
  system.E.coeffRef(this->global_eqn_ids[12], this->global_var_ids[14]) = 1.0;
  // DOF 15, Eq 13: Left ventricle outflow
  system.E.coeffRef(this->global_eqn_ids[13], this->global_var_ids[15]) =
      parameters[this->global_param_ids[ParamId::LLV_A]];
}

template <typename T>
void ClosedLoopHeartPulmonary<T>::update_time(ALGEBRA::SparseSystem<T> &system,
                                              std::vector<T> &parameters) {
  this->get_activation_and_elastance_functions(parameters);
}

template <typename T>
void ClosedLoopHeartPulmonary<T>::update_solution(
    ALGEBRA::SparseSystem<T> &system, std::vector<T> &parameters,
    Eigen::Matrix<T, Eigen::Dynamic, 1> &y,
    Eigen::Matrix<T, Eigen::Dynamic, 1> &dy) {
  this->get_psi_ra_la(parameters, y);
  this->get_valve_positions(y);

  // F and C matrices depend on time and solution
  // Specifying all terms here, including constant terms (which can instead be
  // specified in update_constant) for readability (Doesn't seem to make a
  // difference to compute time) DOF IDs are arranged as inflow
  // [P_in,Q_in,P_out,Q_out,internal variables...]

  // DOF 0, Eq 0: Right atrium pressure
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[0]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[4]) =
      -this->AA * parameters[this->global_param_ids[ParamId::EMAX_RA]];
  system.C(this->global_eqn_ids[0]) =
      this->AA * parameters[this->global_param_ids[ParamId::EMAX_RA]] *
          parameters[this->global_param_ids[ParamId::VASO_RA]] +
      psi_ra * (this->AA - 1.0);
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
      parameters[this->global_param_ids[ParamId::RRA_V]] * valves[5];
  system.F.coeffRef(this->global_eqn_ids[3], this->global_var_ids[0]) = -1.0;
  system.F.coeffRef(this->global_eqn_ids[3], this->global_var_ids[6]) = 1.0;

  // DOF 6, Eq 4: Right ventricle pressure
  system.F.coeffRef(this->global_eqn_ids[4], this->global_var_ids[6]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[4], this->global_var_ids[7]) =
      -this->Erv;
  system.C(this->global_eqn_ids[4]) =
      this->Erv * parameters[this->global_param_ids[ParamId::VRV_U]];

  // DOF 7, Eq 5: Right ventricle volume
  system.F.coeffRef(this->global_eqn_ids[5], this->global_var_ids[5]) =
      -1.0 * valves[5];
  system.F.coeffRef(this->global_eqn_ids[5], this->global_var_ids[8]) =
      1.0 * valves[8];

  // DOF 8, Eq 6: Right ventricle outflow
  system.F.coeffRef(this->global_eqn_ids[6], this->global_var_ids[6]) = -1.0;
  system.F.coeffRef(this->global_eqn_ids[6], this->global_var_ids[9]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[6], this->global_var_ids[8]) =
      parameters[this->global_param_ids[ParamId::RRV_A]] * valves[8];

  // DOF 9, Eq 7: Pulmonary pressure
  system.F.coeffRef(this->global_eqn_ids[7], this->global_var_ids[8]) =
      -valves[8];
  system.F.coeffRef(this->global_eqn_ids[7], this->global_var_ids[9]) =
      1.0 / parameters[this->global_param_ids[ParamId::RPD]];
  system.F.coeffRef(this->global_eqn_ids[7], this->global_var_ids[10]) =
      -1.0 / parameters[this->global_param_ids[ParamId::RPD]];

  // DOF 10, Eq 8: Left atrium pressure
  system.F.coeffRef(this->global_eqn_ids[8], this->global_var_ids[10]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[8], this->global_var_ids[11]) =
      -this->AA * parameters[this->global_param_ids[ParamId::EMAX_LA]];
  system.C(this->global_eqn_ids[8]) =
      this->AA * parameters[this->global_param_ids[ParamId::EMAX_LA]] *
          parameters[this->global_param_ids[ParamId::VASO_LA]] +
      psi_la * (this->AA - 1.0);
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
      parameters[this->global_param_ids[ParamId::RLA_V]] * valves[12];

  // DOF 13, Eq 11: Left ventricle pressure
  system.F.coeffRef(this->global_eqn_ids[11], this->global_var_ids[13]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[11], this->global_var_ids[14]) =
      -this->Elv;
  system.C(this->global_eqn_ids[11]) =
      this->Elv * parameters[this->global_param_ids[ParamId::VLV_U]];

  // DOF 14, Eq 12: Left ventricle volume
  system.F.coeffRef(this->global_eqn_ids[12], this->global_var_ids[12]) =
      -1.0 * valves[12];
  system.F.coeffRef(this->global_eqn_ids[12], this->global_var_ids[15]) =
      1.0 * valves[15];

  // DOF 15, Eq 13: Left ventricle outflow
  system.F.coeffRef(this->global_eqn_ids[13], this->global_var_ids[13]) = -1.0;
  system.F.coeffRef(this->global_eqn_ids[13], this->global_var_ids[2]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[13], this->global_var_ids[15]) =
      parameters[this->global_param_ids[ParamId::RLV_AO]] * valves[15];
}

template <typename T>
void ClosedLoopHeartPulmonary<T>::get_activation_and_elastance_functions(
    std::vector<T> &parameters) {
  T T_cardiac = this->model->cardiac_cycle_period;
  T Tsa = T_cardiac * parameters[this->global_param_ids[ParamId::TSA]];
  T tpwave = T_cardiac / parameters[this->global_param_ids[ParamId::TPWAVE]];
  T t_in_cycle = fmod(this->model->time, T_cardiac);

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

  Elv = Elv_i * parameters[this->global_param_ids[ParamId::ELV_S]];
  Erv = Elv_i * parameters[this->global_param_ids[ParamId::ERV_S]];
}

template <typename T>
void ClosedLoopHeartPulmonary<T>::get_psi_ra_la(
    std::vector<T> &parameters, Eigen::Matrix<T, Eigen::Dynamic, 1> &y) {
  auto RA_volume = y[this->global_var_ids[4]];
  auto LA_volume = y[this->global_var_ids[11]];
  psi_ra =
      parameters[this->global_param_ids[ParamId::KXP_RA]] *
      (exp((RA_volume - parameters[this->global_param_ids[ParamId::VASO_RA]]) *
           parameters[this->global_param_ids[ParamId::KXV_RA]]) -
       1.0);
  psi_la =
      parameters[this->global_param_ids[ParamId::KXP_LA]] *
      (exp((LA_volume - parameters[this->global_param_ids[ParamId::VASO_LA]]) *
           parameters[this->global_param_ids[ParamId::KXV_LA]]) -
       1.0);

  psi_ra_derivative =
      parameters[this->global_param_ids[ParamId::KXP_RA]] *
      exp((RA_volume - parameters[this->global_param_ids[ParamId::VASO_RA]]) *
          parameters[this->global_param_ids[ParamId::KXV_RA]]) *
      parameters[this->global_param_ids[ParamId::KXV_RA]];
  psi_la_derivative =
      parameters[this->global_param_ids[ParamId::KXP_LA]] *
      exp((LA_volume - parameters[this->global_param_ids[ParamId::VASO_LA]]) *
          parameters[this->global_param_ids[ParamId::KXV_LA]]) *
      parameters[this->global_param_ids[ParamId::KXV_LA]];
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
std::map<std::string, int> ClosedLoopHeartPulmonary<T>::get_num_triplets() {
  return num_triplets;
}

}  // namespace MODEL

#endif  // SVZERODSOLVER_MODEL_CLOSEDLOOPHEARTPULMONARY_HPP_
