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
 * ### Parameters
 *
 * * `0` Emax: Maximum elastance
 * * `1` Emin: Minimum elastance
 * * `2` Vrd: Rest diastolic volume
 * * `3` Vrs: Rest systolic volume
 * * `4` Impedance: Impedance of the outflow
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

  void setup_dofs(DOFHandler& dofhandler);
  void update_constant(SparseSystem& system, std::vector<double>& parameters);
  void update_time(SparseSystem& system, std::vector<double>& parameters);

  TripletsContributions num_triplets{6, 2, 0};

  void set_activation_function(std::unique_ptr<ActivationFunction> af) override;

 protected:
  /**
   * @brief Construct a ChamberElastanceInductor with custom block type and
   * parameters. Used by derived classes.
   */
  ChamberElastanceInductor(
      int id, Model* model, BlockType block_type,
      std::vector<std::pair<std::string, InputParameter>> params)
      : Block(id, model, block_type, BlockClass::chamber, params) {}

  double Elas = 0.0;   ///< Current chamber elastance
  double Vrest = 0.0;  ///< Current rest volume
  double act_ = 0.0;   ///< Last computed activation
  std::unique_ptr<ActivationFunction> activation_func_;

  /**
   * @brief Compute elastance and rest volume from activation and parameters.
   */
  virtual void get_elastance_values(std::vector<double>& parameters);
};

#endif  // SVZERODSOLVER_MODEL_CHAMBERELASTANCEINDUCTOR_HPP_
