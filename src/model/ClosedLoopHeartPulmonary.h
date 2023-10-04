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
 * @file ClosedLoopHeartPulmonary.h
 * @brief model::ClosedLoopHeartPulmonary source file
 */
#ifndef SVZERODSOLVER_MODEL_CLOSEDLOOPHEARTPULMONARY_HPP_
#define SVZERODSOLVER_MODEL_CLOSEDLOOPHEARTPULMONARY_HPP_

#include "Block.h"
#include "SparseSystem.h"

// [TODO] get rid of PI.
#define PI 3.14159265
#include <math.h>

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
 */
class ClosedLoopHeartPulmonary : public Block {
 public:
  // Inherit constructors
  using Block::Block;

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
  void update_constant(SparseSystem &system,
                       std::vector<double> &parameters);

  /**
   * @brief Update the time-dependent contributions of the element in a sparse
   * system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  void update_time(SparseSystem &system,
                   std::vector<double> &parameters);

  /**
   * @brief Update the solution-dependent contributions of the element in a
   * sparse system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   * @param y Current solution
   * @param dy Current derivate of the solution
   */
  void update_solution(SparseSystem &system,
                       std::vector<double> &parameters,
                       Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
                       Eigen::Matrix<double, Eigen::Dynamic, 1> &dy);

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
                     Eigen::Matrix<double, Eigen::Dynamic, 1> &y);

  /**
   * @brief Valve positions for each heart chamber
   *
   * @param y Current solution
   */
  void get_valve_positions(Eigen::Matrix<double, Eigen::Dynamic, 1> &y);
};

#endif  // SVZERODSOLVER_MODEL_CLOSEDLOOPHEARTPULMONARY_HPP_
