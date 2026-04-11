// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file BloodVesselFC.h
 * @brief model::BloodVesselFC source file - Blood vessel with fixed capacitance
 */
#ifndef SVZERODSOLVER_MODEL_BLOODVESSELFC_HPP_
#define SVZERODSOLVER_MODEL_BLOODVESSELFC_HPP_

#include <math.h>

#include "Block.h"
#include "SparseSystem.h"

/**
 * @brief Resistor-capacitor-inductor blood vessel with fixed capacitance
 *
 * Same as BloodVessel, but capacitance is not an optimization parameter.
 * Instead, capacitance is read from the model's fixed_capacitance map.
 * This is used in the calibrator when capacitance should not be optimized.
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block (note: no capacitance)
 *
 * * `0` Poiseuille resistance
 * * `1` Inductance
 * * `2` Stenosis coefficient
 *
 */
class BloodVesselFC : public Block {
 public:
  /**
   * @brief Local IDs of the parameters
   *
   */
  enum ParamId {
    RESISTANCE = 0,
    INDUCTANCE = 1,
    STENOSIS_COEFFICIENT = 2,
  };

  /**
   * @brief Construct a new BloodVesselFC object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  BloodVesselFC(int id, Model* model)
      : Block(id, model, BlockType::blood_vessel, BlockClass::vessel,
              {{"R_poiseuille", InputParameter()},
               {"L", InputParameter(true)},
               {"stenosis_coefficient", InputParameter(true)}}) {}

  /**
   * @brief Set up the degrees of freedom (DOF) of the block
   *
   * @param dofhandler Degree-of-freedom handler to register variables and
   * equations at
   */
  void setup_dofs(DOFHandler& dofhandler);

  /**
   * @brief Update the constant contributions of the element in a sparse system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  void update_constant(SparseSystem& system, std::vector<double>& parameters);

  /**
   * @brief Update the solution-dependent contributions of the element in a
   * sparse system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   * @param y Current solution
   * @param dy Current derivate of the solution
   */
  void update_solution(SparseSystem& system, std::vector<double>& parameters,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy);

  /**
   * @brief Set the gradient of the block contributions with respect to the
   * parameters
   *
   * @param jacobian Jacobian with respect to the parameters
   * @param alpha Current parameter vector
   * @param residual Residual with respect to the parameters
   * @param y Current solution
   * @param dy Time-derivative of the current solution
   */
  void update_gradient(Eigen::SparseMatrix<double>& jacobian,
                       Eigen::Matrix<double, Eigen::Dynamic, 1>& residual,
                       Eigen::Matrix<double, Eigen::Dynamic, 1>& alpha,
                       std::vector<double>& y, std::vector<double>& dy);

  /**
   * @brief Number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  TripletsContributions num_triplets{5, 3, 2};
};

#endif  // SVZERODSOLVER_MODEL_BLOODVESSELFC_HPP_
