// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file OpenLoopCoronaryVarResBC.h
 * @brief model::OpenLoopCoronaryVarResBC source file
 */
#ifndef SVZERODSOLVER_MODEL_OPENLOOPCORONARYVARRESBC_HPP_
#define SVZERODSOLVER_MODEL_OPENLOOPCORONARYVARRESBC_HPP_

#include "Block.h"
#include "Parameter.h"
#include "SparseSystem.h"

/**
 * @brief Open loop coronary boundary condition with time-varying microvascular
 * resistance.
 *
 * This model extends the standard open loop coronary BC by allowing the
 * microvascular resistance \f$R_{am}\f$ to vary with time according to:
 *
 * \f[
 * R_{am}(t) = \left[ \left( \sqrt{R_{am,max}} - \sqrt{R_{am,min}} \right)
 * e(t) + \sqrt{R_{am,min}} \right]^2
 * \f]
 *
 * where \f$e(t)\f$ is a time-dependent function:
 *
 * \f[
 * e(t) =
 * \begin{cases}
 * \dfrac{1}{2}\left[ 1 - \cos\left( \pi \dfrac{t_{cycle}}{T_{vc}} \right)
 * \right], & 0 \le t_{cycle} \le T_{vc}, \\
 * \dfrac{1}{2}\left[ 1 + \cos\left( \pi \dfrac{t_{cycle} - T_{vc}}{T_{vr}}
 * \right) \right], & T_{vc} < t_{cycle} \le T_{vc} + T_{vr}, \\
 * 0, & t_{cycle} > T_{vc} + T_{vr}
 * \end{cases}
 * \f]
 *
 * where \f$t_{cycle}\f$ is the time within the current cardiac cycle.
 *
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0) node[anchor=south]{$P_{in}$}
 * to [R, l=$R_a$, *-*] (3,0)
 * to [R, l=$R_{am}(t)$, -] (5,0)
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
 * Same as OpenLoopCoronaryBC, but with \f$R_{am}\f$ replaced by
 * \f$R_{am}(t)\f$.
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block
 *
 * * `0` Ra: Small artery resistance
 * * `1` Ram_min: Minimum microvascular resistance
 * * `2` Ram_max: Maximum microvascular resistance
 * * `3` Rv: Venous resistance
 * * `4` Ca: Small artery capacitance
 * * `5` Cim: Intramyocardial capacitance
 * * `6` Pim: Intramyocardial pressure (time-dependent)
 * * `7` Pv: Venous pressure
 * * `8` T_vc: Contraction time
 * * `9` T_vr: Relaxation time
 *
 * ### Usage in json configuration file
 *
 *     "boundary_conditions": [
 *         {
 *             "bc_name": "OUT",
 *             "bc_type": "CORONARY_VAR_RES",
 *             "bc_values": {
 *                 "Ca": 0.0001,
 *                 "Cc": 0.0001,
 *                 "Pim": [
 *                     1000.0,
 *                     1000.0
 *                 ],
 *                 "P_v": 0.0,
 *                 "Ra1": 100.0,
 *                 "Ra2_min": 50.0,
 *                 "Ra2_max": 200.0,
 *                 "Rv1": 100.0,
 *                 "T_vc": 0.2,
 *                 "T_vr": 0.3,
 *                 "t": [
 *                     0.0,
 *                     1.0
 *                 ]
 *             }
 *         }
 *     ]
 *
 * ### Internal variables
 *
 * Names of internal variables in this block's output:
 *
 * * `volume_im`: Intramyocardial volume
 *
 */
class OpenLoopCoronaryVarResBC : public Block {
 public:
  /**
   * @brief Construct a new OpenLoopCoronaryVarResBC object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  OpenLoopCoronaryVarResBC(int id, Model* model)
      : Block(id, model, BlockType::open_loop_coronary_var_res_bc,
              BlockClass::boundary_condition,
              {{"Ra1", InputParameter()},
               {"Ra2_min", InputParameter()},
               {"Ra2_max", InputParameter()},
               {"Rv1", InputParameter()},
               {"Ca", InputParameter()},
               {"Cc", InputParameter()},
               {"t", InputParameter(false, true)},
               {"Pim", InputParameter(false, true)},
               {"P_v", InputParameter()},
               {"T_vc", InputParameter()},
               {"T_vr", InputParameter()},
               {"closed_loop_outlet", InputParameter(true, false, false)}}) {}

  /**
   * @brief Set up the degrees of freedom (DOF) of the block
   *
   * Set \ref global_var_ids and \ref global_eqn_ids of the element based on
   * the number of equations and the number of internal variables of the
   * element.
   *
   * @param dofhandler Degree-of-freedom handler to register variables and
   * equations at
   */
  void setup_dofs(DOFHandler& dofhandler);

  /**
   * @brief Setup parameters that depend on the initial state
   *
   * @param initial_state The initial state of the system
   * @param parameters The parameter values vector (at time 0)
   */
  void setup_initial_state_dependent_params(State initial_state,
                                            std::vector<double>& parameters);

  /**
   * @brief Update the constant contributions of the element in a sparse system
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
  TripletsContributions num_triplets{5, 4, 0};

 protected:
  double P_Cim_0 = 0;  ///< Pressure proximal to Cim/Vim at initial state
  double Pim_0 = 0;    ///< Pim at initial state

  /**
   * @brief Compute time-varying microvascular resistance
   *
   * @param t_cycle Current time within cardiac cycle
   * @param Ram_min Minimum microvascular resistance
   * @param Ram_max Maximum microvascular resistance
   * @param T_vc Contraction time
   * @param T_vr Relaxation time
   * @return Time-varying resistance Ram(t)
   */
  double compute_Ram(double t_cycle, double Ram_min, double Ram_max,
                     double T_vc, double T_vr) const;
};

#endif  // SVZERODSOLVER_MODEL_OPENLOOPCORONARYVARRESBC_HPP_
