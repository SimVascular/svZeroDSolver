// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

/**
 * @file LinearElastanceChamber.h
 * @brief model::LinearElastanceChamber source file
 */

#ifndef SVZERODSOLVER_MODEL_LINEARELASTANCECHAMBER_HPP_
#define SVZERODSOLVER_MODEL_LINEARELASTANCECHAMBER_HPP_

#include <math.h>

#include <memory>

#include "ActivationFunction.h"
#include "Block.h"
#include "Model.h"
#include "SparseSystem.h"
#include "debug.h"

/**
 * @brief Cardiac chamber with linear elastance (no inductor).
 *
 * Models a cardiac chamber as a time-varying capacitor (elastance with
 * specified resting volumes) without inductance. See \cite Regazzoni2022
 *
 * This chamber block can be connected to other blocks using junctions.
 *
 * \f[
 * \begin{circuitikz}
 * \draw node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0) node[anchor=south]{$P_{in}$}
 * to[short, *-*] (3,0) node[anchor=south]{$P_{out}$}
 * (1,0) to [vC, l=$E$, *-] (1,-1.5)
 * node[ground]{};
 * \draw (3.2,0) -- (4.0,0) node[right] {$Q_{out}$} [-latex];
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f[
 * P_{in}-E(t)(V_c-V_{rest})=0
 * \f]
 *
 * \f[
 * P_{in}-P_{out}=0
 * \f]
 *
 * \f[
 * Q_{in}-Q_{out}-\dot{V}_c=0
 * \f]
 *
 * ### Local contributions
 *
 * \f[
 * \mathbf{y}^{e}=\left[\begin{array}{lllll}P_{in} & Q_{in} &
 * P_{out} & Q_{out} & V_c\end{array}\right]^{T} \f]
 *
 * \f[
 * \mathbf{E}^{e}=\left[\begin{array}{ccccc}
 * 0 & 0 & 0 & 0 & 0\\
 * 0 & 0 & 0 & 0 & 0\\
 * 0 & 0 & 0 & 0 & -1
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{F}^{e}=\left[\begin{array}{ccccc}
 * 1 & 0 &  0 & 0  & E(t) \\
 * 1 & 0 & -1 & 0  & 0 \\
 * 0 & 1 &  0 & -1 & 0
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{c}^{e}=\left[\begin{array}{c}
 * E(t)V_{rest} \\
 * 0 \\
 * 0
 * \end{array}\right]
 * \f]
 *
 * In the above equations,
 *
 * \f[
 * E_i(t) = E_i^{\text{pass}} + E_i^{\text{act,max}} \,
 * \phi\!\left(t, t_C^i, t_R^i, T_C^i, T_R^i\right),
 * \f]
 *
 * \f[
 * \phi(t, t_C, t_R, T_C, T_R) =
 * \begin{cases}
 * \frac{1}{2}\left[1 - \cos\!\left(\frac{\pi}{T_C} \operatorname{mod}(t - t_C,
 * T_{\mathrm{HB}})\right)\right], & \text{if } 0 \le \operatorname{mod}(t -
 * t_C, T_{\mathrm{HB}}) < T_C, \\[1.2em] \frac{1}{2}\left[1 +
 * \cos\!\left(\frac{\pi}{T_R} \operatorname{mod}(t - t_R,
 * T_{\mathrm{HB}})\right)\right], & \text{if } 0 \le \operatorname{mod}(t -
 * t_R, T_{\mathrm{HB}}) < T_R, \\[1.2em] 0, & \text{otherwise.} \end{cases} \f]
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block
 *
 * * `0` Emax: Maximum elastance
 * * `1` Epass: Passive elastance
 * * `2` Vrest: Rest diastolic volume
 * * `3` contract_start: Contract start time
 * * `4` relax_start: Relax start time
 * * `5` contract_duration: Contract duration
 * * `6` relax_duration: Relaxation duration
 *
 * ### Internal variables
 *
 * Names of internal variables in this block's output:
 *
 * * `Vc`: Chamber volume
 *
 */
class LinearElastanceChamber : public Block {
 public:
  /**
   * @brief Construct a new LinearElastanceChamber object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  LinearElastanceChamber(int id, Model* model)
      : Block(id, model, BlockType::linear_elastance_chamber,
              BlockClass::chamber,
              {{"Emax", InputParameter()},
               {"Epass", InputParameter()},
               {"Vrest", InputParameter()}}) {}

  /**
   * @brief Local IDs of the parameters
   *
   */
  enum ParamId { EMAX = 0, EPASS = 1, VREST = 2 };

  /**
   * @brief Set up the degrees of freedom (DOF) of the block
   *
   * Set global_var_ids and global_eqn_ids of the element based on the
   * number of equations and the number of internal variables of the
   * element.
   *
   * @param dofhandler Degree-of-freedom handler to register variables and
   * equations at
   */
  void setup_dofs(DOFHandler& dofhandler);

  /**
   * @brief Update the constant contributions of the element in a sparse
   system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  void update_constant(SparseSystem& system, std::vector<double>& parameters);

  /**
   * @brief Update the time-dependent contributions of the element in a sparse
   * system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  void update_time(SparseSystem& system, std::vector<double>& parameters);

  /**
   * @brief Number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  TripletsContributions num_triplets{6, 2, 0};

 private:
  double Elas;                                           // Chamber Elastance
  std::unique_ptr<ActivationFunction> activation_func_;  // Activation function

 public:
  /**
   * @brief Set the activation function (takes ownership)
   *
   * @param af Unique pointer to the activation function
   */
  void set_activation_function(std::unique_ptr<ActivationFunction> af) override;

 private:
  /**
   * @brief Update the elastance functions which depend on time
   *
   * @param parameters Parameters of the model
   */
  void get_elastance_values(std::vector<double>& parameters);
};

#endif  // SVZERODSOLVER_MODEL_LINEARELASTANCECHAMBER_HPP_
