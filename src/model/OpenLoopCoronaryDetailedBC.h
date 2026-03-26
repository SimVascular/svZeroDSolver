// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file OpenLoopCoronaryDetailedBC.h
 * @brief model::OpenLoopCoronaryDetailedBC source file
 */
#ifndef SVZERODSOLVER_MODEL_OPENLOOPCORONARYDETAILEDBC_HPP_
#define SVZERODSOLVER_MODEL_OPENLOOPCORONARYDETAILEDBC_HPP_

#include "Block.h"
#include "Parameter.h"
#include "SparseSystem.h"

/**
 * @brief Open loop coronary boundary condition with detailed internal outputs.
 *
 * This model extends the standard open loop coronary BC to output internal
 * pressures and flows at each segment of the coronary circuit.
 *
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0) node[anchor=south]{$P_{in}$}
 * to [R, l=$R_{a1}$, *-*] (3,0)
 * node[anchor=south]{$P_a$}
 * to [R, l=$R_{a2}$, -] (5,0)
 * node[anchor=south]{$P_{c}$}
 * to [R, l=$R_{v1}$, *-*] (7,0)
 * node[anchor=south]{$P_{v}$}
 * (5,0) to [C, l=$C_c \;V_c$, -*] (5,-1.5)
 * node[left]{$P_{im}$}
 * (3,0) to [C, l=$C_a \; V_a$, -*] (3,-1.5)
 * node[ground]{};
 * \end{circuitikz}
 * \f]
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block
 *
 * * `0` Ra1: Arterial resistance
 * * `1` Rv1: Venous resistance
 * * `2` Ca: Arterial capacitance
 * * `3` Cc: Capillary capacitance
 * * `4` Pim: Intramyocardial pressure (time-dependent)
 * * `5` Pv: Venous pressure
 * * `6` Ra2: Microvascular resistance
 *
 * ### Usage in json configuration file
 *
 *     "boundary_conditions": [
 *         {
 *             "bc_name": "OUT",
 *             "bc_type": "CORONARY_DETAILED",
 *             "bc_values": {
 *                 "Ca": 0.0001,
 *                 "Cc": 0.0001,
 *                 "Pim": [1000.0, 1000.0],
 *                 "P_v": 0.0,
 *                 "Ra1": 100.0,
 *                 "Ra2": 100.0,
 *                 "Rv1": 100.0,
 *                 "t": [0.0, 1.0]
 *             }
 *         }
 *     ]
 *
 * ### Internal variables
 *
 * Names of internal variables in this block's output:
 *
 * * `V_a`: Arterial compliance volume (Ca * Pa)
 * * `V_c`: Capillary compliance volume (Cc * (Pc - Pim))
 * * `P_a`: Pressure after Ra1 (at arterial compliance)
 * * `P_c`: Pressure after Ra2 (at capillary compliance)
 *
 */
class OpenLoopCoronaryDetailedBC : public Block {
 public:
  /**
   * @brief Construct a new OpenLoopCoronaryDetailedBC object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  OpenLoopCoronaryDetailedBC(int id, Model* model)
      : Block(id, model, BlockType::open_loop_coronary_detailed_bc,
              BlockClass::boundary_condition,
              {{"Ra1", InputParameter()},
               {"Rv1", InputParameter()},
               {"Ca", InputParameter()},
               {"Cc", InputParameter()},
               {"t", InputParameter(false, true)},
               {"Pim", InputParameter(false, true)},
               {"P_v", InputParameter()},
               {"Ra2", InputParameter()},
               {"closed_loop_outlet", InputParameter(true, false, false)}}) {}

  /**
   * @brief Set up the degrees of freedom (DOF) of the block
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
   * E matrix: 2 entries (dV_a/dt, dV_c/dt)
   * F matrix: 15 entries across 5 equations
   */
  TripletsContributions num_triplets{15, 2, 0};

 protected:
  double Pim_0 = 0;  ///< Pim at initial state
};

#endif  // SVZERODSOLVER_MODEL_OPENLOOPCORONARYDETAILEDBC_HPP_
