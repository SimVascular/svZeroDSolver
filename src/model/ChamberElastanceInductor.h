// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file ChamberElastanceInductor.h
 * @brief model::ChamberElastanceInductor source file
 */
#ifndef SVZERODSOLVER_MODEL_CHAMBERELASTANCEINDUCTOR_HPP_
#define SVZERODSOLVER_MODEL_CHAMBERELASTANCEINDUCTOR_HPP_

#include <math.h>

#include <memory>

#include "ActivationFunction.h"
#include "Block.h"
#include "SparseSystem.h"
#include "debug.h"

/**
 * @brief Cardiac chamber with linear elastance and inductor.
 *
 * Models a cardiac chamber as a time-varying capacitor (elastance with
 * specified resting volumes) and an inductor. See \cite kerckhoffs2007coupling
 * (equations 1 and 2). The addition of the inductor is similar to the models in
 * \cite sankaran2012patient and \cite menon2023predictors.
 *
 * This chamber block can be connected to other blocks using junctions.
 *
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0) node[anchor=south]{$P_{in}$}
 * to (1,0)
 * node[anchor=south]{}
 * to [L, l=$L$, *-*] (3,0)
 * node[anchor=south]{$P_{out}$}
 * (1,0) to [vC, l=$E$, *-] (1,-1.5)
 * node[ground]{};
 * \draw [-latex] (3.2,0) -- (4.0,0) node[right] {$Q_{out}$} ;
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
 * P_{in}-P_{out}-L\dot{Q}_{out}=0
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
 * 0 & 0 & 0 & -L & 0\\
 * 0 & 0 & 0 & 0 & -1
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{F}^{e}=\left[\begin{array}{ccccc}
 * 1 & 0 &  0 & 0  & -E(t) \\
 * 1 & 0 & -1 & 0  & 0 \\
 * 0 & 1 &  0 & -1 & 0
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{c}^{e}=\left[\begin{array}{c}
 * E(t) V_{rest} \\
 * 0 \\
 * 0
 * \end{array}\right]
 * \f]
 *
 * where
 *
 * \f[
 * V_{rest}(t)= \{1-A(t)\}(V_{rd}-V_{rs})+V_{rs}
 * \f]
 *
 * \f[
 * E(t)=(E_{max}-E_{min})A(t) + E_{min}
 * \f]
 *
 * The activation function \f$A(t)\f$ is provided as a separate object
 * (e.g. half_cosine, fourier, wrapping_cosine).
 *
 * ### Parameters
 *
 * * `0` Impedance: Outflow inductance \f$L\f$
 * * `1` Emax: Maximum elastance
 * * `2` Emin: Minimum elastance
 * * `3` Vrd: Rest diastolic volume
 * * `4` Vrs: Rest systolic volume
 *
 * ### Internal variables
 *
 * * `Vc`: Chamber volume
 *
 */
class ChamberElastanceInductor : public Block {
 public:
  /**
   * @brief Construct a new ChamberElastanceInductor object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  ChamberElastanceInductor(int id, Model* model)
      : Block(id, model, BlockType::chamber_elastance_inductor,
              BlockClass::chamber,
              {{"Impedance", InputParameter()},
               {"Emax", InputParameter()},
               {"Emin", InputParameter()},
               {"Vrd", InputParameter()},
               {"Vrs", InputParameter()}}) {}

  /**
   * @brief Local IDs of the parameters (shared indices first)
   */
  enum ParamId {
    IMPEDANCE = 0,
    EMAX = 1,
    EMIN = 2,
    VRD = 3,
    VRS = 4,
  };

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
  void setup_dofs(DOFHandler& dofhandler) override;

  /**
   * @brief Update the constant contributions of the element in a sparse
   system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  void update_constant(SparseSystem& system,
                       std::vector<double>& parameters) override;

  /**
   * @brief Update the time-dependent contributions of the element in a sparse
   * system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  void update_time(SparseSystem& system,
                   std::vector<double>& parameters) override;

  /// @brief Number of triplets of element
  TripletsContributions num_triplets{6, 2, 0};

  void set_activation_function(std::unique_ptr<ActivationFunction> af) override;

 protected:
  /**
   * @brief Construct a ChamberElastanceInductor with custom block type and
   * parameters. Used by derived classes.
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   * @param block_type The specific block type for the derived class
   * @param params Parameter list for the derived class
   */
  ChamberElastanceInductor(
      int id, Model* model, BlockType block_type,
      std::vector<std::pair<std::string, InputParameter>> params)
      : Block(id, model, block_type, BlockClass::chamber, params) {}

  double Elas = 0.0;   ///< Current chamber elastance
  double Vrest = 0.0;  ///< Current rest volume
  double act_ = 0.0;   ///< Last computed activation
  /// @brief Activation function
  std::unique_ptr<ActivationFunction> activation_func_;

  /**
   * @brief Compute elastance and rest volume from activation and parameters.
   *
   * @param parameters Parameters of the model
   */
  virtual void get_elastance_values(std::vector<double>& parameters);
};

#endif  // SVZERODSOLVER_MODEL_CHAMBERELASTANCEINDUCTOR_HPP_
