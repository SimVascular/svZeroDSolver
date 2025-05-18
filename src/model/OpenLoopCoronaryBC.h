// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file OpenLoopCoronaryBC.h
 * @brief model::OpenLoopCoronaryBC source file
 */
#ifndef SVZERODSOLVER_MODEL_OPENLOOPCORONARYBC_HPP_
#define SVZERODSOLVER_MODEL_OPENLOOPCORONARYBC_HPP_

#include "Block.h"
#include "Parameter.h"
#include "SparseSystem.h"

/**
 * @brief Open loop coronary boundary condition based on \cite kim_coronary.
 *
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0) node[anchor=south]{$P_{in}$}
 * to [R, l=$R_a$, *-*] (3,0)
 * to [R, l=$R_{am}$, -] (5,0)
 * node[anchor=south]{$P_{cim}$}
 * to [R, l=$R_v$, *-*] (7,0)
 * node[anchor=south]{$P_{v}$}
 * (5,0) to [C, l=$C_{im} \;V_{im}$, -*] (5,-1.5)
 * node[left]{$P_{im}$}
 * (3,0) to [C, l=$C_a$, -*] (3,-1.5)
 * node[left]{$P_a$};
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f{eqnarray*}{
 * &C_{i m} R_{v} Q_{in}-V_{i m}+C_{i m} \left(-P_{c i m}(0)+P_{i m}(0)-P_{i
 * m}+P_{v}\right)-C_{i m} R_{v} \frac{d V_{i m}}{d t}-C_{a} C_{i m} R_{v}
 * \frac{d P_{in}}{d t}+R_{a} C_{a} C_{i m} R_{v} \frac{d Q_{in}}{d t}\\ &+C_{a}
 * C_{i m} R_{v} \frac{d P_{a}}{d t}=0
 * \f}
 *
 * \f{eqnarray*}{
 * &C_{i m} R_v P_{in}-C_{i m} R_{v} R_{a} Q_{in}-R_{v} V_{i m}-C_{i
 * m}\left(R_{v}+R_{a m}\right) \left(P_{c i m}(0)-P_{i m}(0)+P_{i
 * m}\right)-C_{i m} R_{v} R_{a m} \frac{d V_{i m}}{d t} \\ &-R_{a m} V_{i
 * m}+R_{a m} C_{i m} P_{v}=0 \f}
 *
 * ### Local contributions
 *
 * \f[
 * \mathbf{y}^{e}=\left[\begin{array}{lll}P_{in} & Q_{in} & V_{i
 * m}\end{array}\right]^{T}, \f]
 *
 * \f[
 * \mathbf{E}^{e}=\left[\begin{array}{ccc}-C_{a} C_{i m} R_{v} & R_{a} C_{a}
 * C_{i m} R_{v} & -C_{i m} R_{v} \\ 0 & 0 & -C_{i m} R_{v} R_{a
 * m}\end{array}\right] \f]
 *
 * \f[
 * \mathbf{F}^{e}=\left[\begin{array}{ccc}0 & C_{i m} R_{v} & -1 \\C_{i m} R_{v}
 * & -C_{i m} R_{v} R_{a} & -\left(R_{v}+R_{a m}\right)\end{array}\right] \f]
 *
 * \f[
 * \mathbf{c}^{e}=\left[\begin{array}{c}C_{i m}\left(-P_{i m}+P_{i m}(0)-P_{c i
 * m}(0)+P_{v}\right)+C_{a} C_{i m} R_{v} \frac{d P_{a}}{d t} \\-C_{i
 * m}\left(R_{v} + R_{a m}\right)\left(P_{cim}(0)-P_{i m}(0)+P_{i m}\right)+R_{a
 * m} C_{i m} P_{v}\end{array}\right] \f]
 *
 * Assume \f$P_a=0\f$.
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block
 *
 * * `0` Ra: Small artery resistance
 * * `1` Ram: Microvascualar resistance
 * * `2` Rv: Venous resistance
 * * `3` Ca: Small artery capacitance
 * * `4` Cim: Intramyocardial capacitance
 * * `5` Pim: Intramyocardial pressure
 * * `6` Pv: Venous pressure
 *
 * ### Internal variables
 *
 * Names of internal variables in this block's output:
 *
 * * `volume_im`: Intramyocardial volume
 *
 */
class OpenLoopCoronaryBC : public Block {
 public:
  /**
   * @brief Construct a new OpenLoopCoronaryBC object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  OpenLoopCoronaryBC(int id, Model *model)
      : Block(id, model, BlockType::open_loop_coronary_bc,
              BlockClass::boundary_condition,
              {{"Ra1", InputParameter()},
               {"Ra2", InputParameter()},
               {"Rv1", InputParameter()},
               {"Ca", InputParameter()},
               {"Cc", InputParameter()},
               {"t", InputParameter(false, true)},
               {"Pim", InputParameter(false, true)},
               {"P_v", InputParameter()},
               {"closed_loop_outlet", InputParameter(true, false, false)}}) {}

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
   * @brief Setup parameters that depend on the initial state
   *
   * @param initial_state The initial state of the system
   * @param parameters The parameter values vector (at time 0)
   */
  void setup_initial_state_dependent_params(State initial_state,
                                            std::vector<double> &parameters);

  /**
   * @brief Update the constant contributions of the element in a sparse system
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
   * @brief Number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  TripletsContributions num_triplets{5, 4, 0};

 protected:
  double P_Cim_0 = 0;  ///< Pressure proximal to Cim/Vim at initial state
  double Pim_0 = 0;    ///< Pim at initial state
};

#endif  // SVZERODSOLVER_MODEL_OPENLOOPCORONARYBC_HPP_
