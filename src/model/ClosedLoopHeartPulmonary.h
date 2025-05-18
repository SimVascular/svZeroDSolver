// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file ClosedLoopHeartPulmonary.h
 * @brief model::ClosedLoopHeartPulmonary source file
 */
#ifndef SVZERODSOLVER_MODEL_CLOSEDLOOPHEARTPULMONARY_HPP_
#define SVZERODSOLVER_MODEL_CLOSEDLOOPHEARTPULMONARY_HPP_

#include "Block.h"
#include "SparseSystem.h"

/**
 * @brief Define math constants (for M_PI)
 *
 */
#define _USE_MATH_DEFINES

#include <cmath>

/**
 * @brief Heart and pulmonary circulation model
 *
 * Models the mechanics of the 4 heart chambers and pulmonary circulation
 *
 * References: \cite sankaran2012patient and \cite menon2023predictors
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
 * ### Internal variables
 *
 * Names of internal variables in this block's output:
 *
 * * `V_RA`: Right atrium volume
 * * `Q_RA`: Right atrium outflow
 * * `P_RV`: Right ventricle pressure
 * * `V_RV`: Right ventricle volume
 * * `Q_RV`: Right ventricle outflow
 * * `P_pul`: Pulmonary pressure
 * * `P_LA`: Left atrium pressure
 * * `V_LA`: Left atrium volume
 * * `Q_LA`: Left atrium outflow
 * * `P_LV`: Left ventricle pressure
 * * `V_LV`: Left ventricle volume
 * * `Q_LV`: Left ventricle outflow
 *
 */
class ClosedLoopHeartPulmonary : public Block {
 public:
  /**
   * @brief Construct a new ClosedLoopHeartPulmonary object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  ClosedLoopHeartPulmonary(int id, Model *model)
      : Block(id, model, BlockType::closed_loop_heart_pulmonary,
              BlockClass::closed_loop,
              {{"Tsa", InputParameter()},     {"tpwave", InputParameter()},
               {"Erv_s", InputParameter()},   {"Elv_s", InputParameter()},
               {"iml", InputParameter()},     {"imr", InputParameter()},
               {"Lra_v", InputParameter()},   {"Rra_v", InputParameter()},
               {"Lrv_a", InputParameter()},   {"Rrv_a", InputParameter()},
               {"Lla_v", InputParameter()},   {"Rla_v", InputParameter()},
               {"Llv_a", InputParameter()},   {"Rlv_ao", InputParameter()},
               {"Vrv_u", InputParameter()},   {"Vlv_u", InputParameter()},
               {"Rpd", InputParameter()},     {"Cp", InputParameter()},
               {"Cpa", InputParameter()},     {"Kxp_ra", InputParameter()},
               {"Kxv_ra", InputParameter()},  {"Kxp_la", InputParameter()},
               {"Kxv_la", InputParameter()},  {"Emax_ra", InputParameter()},
               {"Emax_la", InputParameter()}, {"Vaso_ra", InputParameter()},
               {"Vaso_la", InputParameter()}}) {}

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
  void update_constant(SparseSystem &system, std::vector<double> &parameters);

  /**
   * @brief Update the time-dependent contributions of the element in a sparse
   * system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  void update_time(SparseSystem &system, std::vector<double> &parameters);

  /**
   * @brief Update the solution-dependent contributions of the element in a
   * sparse system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   * @param y Current solution
   * @param dy Current derivate of the solution
   */
  void update_solution(SparseSystem &system, std::vector<double> &parameters,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy);

  /**
   * @brief Modify the solution after solving it
   *
   * @param y Current solution
   */
  void post_solve(Eigen::Matrix<double, Eigen::Dynamic, 1> &y);

  /**
   * @brief Number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  TripletsContributions num_triplets{33, 10, 2};

 private:
  // Below variables change every timestep and are then combined with
  // expressions that are updated with solution
  double AA;   // Atrial activation function
  double Elv;  // LV elastance
  double Erv;  // RV elastance
  double psi_ra, psi_la, psi_ra_derivative,
      psi_la_derivative;  // Expressions for atrial activation
  double valves[16];

  /**
   * @brief Update the atrial activation and LV/RV elastance functions which
   * depend on time
   *
   * @param parameters Parameters of the model
   */
  void get_activation_and_elastance_functions(std::vector<double> &parameters);

  /**
   * @brief Compute sub-expressions that are part of atrial elastance and
   * depends on atrial volume from the solution vector
   *
   * @param parameters Parameters of the model
   * @param y Current solution
   */
  void get_psi_ra_la(std::vector<double> &parameters,
                     const Eigen::Matrix<double, Eigen::Dynamic, 1> &y);

  /**
   * @brief Valve positions for each heart chamber
   *
   * @param y Current solution
   */
  void get_valve_positions(const Eigen::Matrix<double, Eigen::Dynamic, 1> &y);
};

#endif  // SVZERODSOLVER_MODEL_CLOSEDLOOPHEARTPULMONARY_HPP_
