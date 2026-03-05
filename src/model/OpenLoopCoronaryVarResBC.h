// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file OpenLoopCoronaryVarResBC.h
 * @brief model::OpenLoopCoronaryVarResBC source file
 */
#ifndef SVZERODSOLVER_MODEL_OPENLOOPCORONARYVARRESBC_HPP_
#define SVZERODSOLVER_MODEL_OPENLOOPCORONARYVARRESBC_HPP_

#include "Block.h"
#include "OpenLoopCoronaryBC.h"
#include "Parameter.h"

/**
 * @brief Open loop coronary boundary condition with time-varying microvascular
 * resistance based on \cite yong25.
 *
 * This model extends the standard open loop coronary BC (\cite kim_coronary)
 * by allowing the microvascular resistance \f$R_{am}\f$ to vary with time
 * according to:
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
 * * `1` Rv: Venous resistance
 * * `2` Ca: Small artery capacitance
 * * `3` Cim: Intramyocardial capacitance
 * * `4` Pim: Intramyocardial pressure (time-dependent)
 * * `5` Pv: Venous pressure
 * * `6` Ram_min: Minimum microvascular resistance
 * * `7` Ram_max: Maximum microvascular resistance
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
class OpenLoopCoronaryVarResBC : public OpenLoopCoronaryBC {
 public:
  /**
   * @brief Construct a new OpenLoopCoronaryVarResBC object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  OpenLoopCoronaryVarResBC(int id, Model* model)
      : OpenLoopCoronaryBC(
            id, model, BlockType::open_loop_coronary_var_res_bc,
            {{"Ra1", InputParameter()},
             {"Rv1", InputParameter()},
             {"Ca", InputParameter()},
             {"Cc", InputParameter()},
             {"t", InputParameter(false, true)},
             {"Pim", InputParameter(false, true)},
             {"P_v", InputParameter()},
             {"Ra2_min", InputParameter()},
             {"Ra2_max", InputParameter()},
             {"T_vc", InputParameter()},
             {"T_vr", InputParameter()},
             {"closed_loop_outlet", InputParameter(true, false, false)}}) {}

 protected:
  /**
   * @brief Get time-varying microvascular resistance
   *
   * Overrides base class to provide time-varying resistance based on
   * cardiac cycle phase.
   *
   * @param parameters Parameters of the model
   * @param time Current simulation time
   * @return Time-varying resistance Ram(t)
   */
  double get_Ram(std::vector<double>& parameters, double time) const override;
};

#endif  // SVZERODSOLVER_MODEL_OPENLOOPCORONARYVARRESBC_HPP_
