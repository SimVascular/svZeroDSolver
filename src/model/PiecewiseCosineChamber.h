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
 * @file PiecewiseCosineChamber.h
 * @brief model::PiecewiseCosineChamber source file
 */
#ifndef SVZERODSOLVER_MODEL_PiecewiseCosineChamber_HPP_
#define SVZERODSOLVER_MODEL_PiecewiseCosineChamber_HPP_

#include <math.h>

#include "Block.h"
#include "Model.h"
#include "SparseSystem.h"
#include "debug.h"

/**
 * @brief Cardiac chamber with time-varying elastance (0D model).
 *
 * Models a cardiac chamber as a time-varying capacitance (elastance) with an
 * unstressed (reference) volume. Pressures are given by an elastance law,
 * and inter-compartment flows are set by valve “diodes” whose resistance
 * switches between R_min (forward flow) and R_max (closed valve), consistent
 * with the closed-loop 0D formulation shown in the figures (Eqs. (6)–(7f), (8)–(9)).
 *
 * This chamber block connects to the rest of the circulation through junctions
 * and valves; the valves obey the piecewise resistance below.
 *
 * ### Diagram (conceptual)
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0) node[anchor=south]{$p_{in}$}
 * to (1,0)
 * node[anchor=south]{}
 * to [generic, l=$R_i(p_{in},p_{out})$, *-*] (3,0)
 * node[anchor=south]{$p_{out}$}
 * (1,0) to [vC, l=$E_i(t)$, *-] (1,-1.5)
 * node[ground]{};
 * \draw [-latex] (3.2,0) -- (4.0,0) node[right] {$Q_{out}$} ;
 * \end{circuitikz}
 * \f]
 *
 * ### Governing relations for a chamber \(i\in\{\text{LA},\text{LV},\text{RA},\text{RV}\}\)
 *
 * **Elastance (pressure–volume law):**
 * \f[
 * p_i(t) \;=\; p_{\mathrm{EX}}(t) \;+\; E_i(t)\,\big(V_i(t) - V_{0,i}\big).
 * \tag{7a–7d form}
 * \f]
 *
 * **Mass conservation in the chamber:**
 * \f[
 * \dot V_i(t) \;=\; Q_{in,i}(t) - Q_{out,i}(t).
 * \f]
 *
 * **Valve (non-ideal diode) flow laws (examples):**
 * \f[
 * Q_{\mathrm{MV}}(t) = \dfrac{p_{\mathrm{LA}}(t) - p_{\mathrm{LV}}(t)}
 * {R_{\mathrm{MV}}\!\big(p_{\mathrm{LA}},p_{\mathrm{LV}}\big)},\quad
 * Q_{\mathrm{AV}}(t) = \dfrac{p_{\mathrm{LV}}(t) - p_{\mathrm{AR}}^{\mathrm{SYS}}(t)}
 * {R_{\mathrm{AV}}\!\big(p_{\mathrm{LV}},p_{\mathrm{AR}}^{\mathrm{SYS}}\big)},
 * \tag{7e–7f form}
 * \f]
 * with analogous definitions for \(Q_{\mathrm{TV}}\) and \(Q_{\mathrm{PV}}\).
 *
 * **Valve resistance switching (for } i\in\{\mathrm{MV},\mathrm{AV},\mathrm{TV},\mathrm{PV}\}):**
 * \f[
 * R_i(p_1,p_2) =
 * \begin{cases}
 * R_{\min}, & p_1 < p_2 \quad \text{(forward/open)} \\
 * R_{\max}, & p_1 \ge p_2 \quad \text{(closed/backward)}
 * \end{cases}
 * \f]
 * where \(p_1\) is the upstream pressure (ahead of the leaflets in the flow
 * direction) and \(p_2\) is downstream. Choose \(R_{\max}\) large so that
 * leakage when closed is negligible; \(R_{\min}>0\) adds small forward losses.
 *
 * ### Chamber activation and time-varying elastance
 *
 * Piecewise activation \(\phi(t,\;t_C,\;t_R,\;T_C,\;T_R)\) over the heartbeat
 * period \(T_{HB}\) (Eq. (8)):
 * \f[
 * \phi(t,t_C,t_R,T_C,T_R) =
 * \begin{cases}
 * \tfrac{1}{2}\Big[1-\cos\!\big(\tfrac{\pi}{T_C}\,\mathrm{mod}(t-t_C,T_{HB})\big)\Big],
 * & 0 \le \mathrm{mod}(t-t_C,T_{HB}) < T_C,\\[6pt]
 * \tfrac{1}{2}\Big[1+\cos\!\big(\tfrac{\pi}{T_R}\,\mathrm{mod}(t-t_R,T_{HB})\big)\Big],
 * & 0 \le \mathrm{mod}(t-t_R,T_{HB}) < T_R,\\[6pt]
 * 0, & \text{otherwise.}
 * \end{cases}
 * \tag{8}
 * \f]
 *
 * Time-varying elastance (Eq. (9)):
 * \f[
 * E_i(t) \;=\; E_i^{\mathrm{pass}} \;+\; E_{i}^{\mathrm{act,max}}\;\phi\!\left(
 * t,\;t_C^i,\;t_R^i,\;T_C^i,\;T_R^i\right).
 * \tag{9}
 * \f]
 *
 * Here \(V_{0,i}\) is the chamber’s unstressed volume; \(p_{\mathrm{EX}}(t)\)
 * is extracardiac/pericardial pressure (often set to \(0\) for simplicity).
 *
 * ### Local unknowns (example ordering for a chamber element)
 * \f[
 * \mathbf{y}^{e}=\big[p_{in}\;\;Q_{in}\;\;p_{out}\;\;Q_{out}\;\;V_i\big]^T.
 * \f]
 * In a minimal (no-inductor) chamber, the residuals are formed from:
 * \f[
 * \begin{aligned}
 * &p_{in}-p_{\mathrm{EX}}-E_i(t)(V_i-V_{0,i})=0,\\
 * &Q_{in}-Q_{out}-\dot V_i=0,\\
 * &Q_{out}-\dfrac{p_{in}-p_{out}}{R_i(p_{in},p_{out})}=0,
 * \end{aligned}
 * \f]
 * which can be assembled into your global system. (If you maintain an
 * inductor model elsewhere, add the corresponding \(L\,\dot Q\) term there.)
 *
 * ### Parameters
 * Provide parameters per chamber \(i\):
 * - `E_pass`  (\(E_i^{\mathrm{pass}}\)) – passive elastance
 * - `E_act_max` (\(E_i^{\mathrm{act,max}}\)) – active elastance amplitude
 * - `V0` (\(V_{0,i}\)) – unstressed volume
 * - `tC`, `tR` – start times of contraction and relaxation
 * - `TC`, `TR` – durations of contraction and relaxation
 * - `T_HB` – heartbeat period
 * - `p_ex` (optional) – extracardiac pressure (default 0)
 *
 * Valve parameters (used by the adjoining valve elements):
 * - `Rmin`, `Rmax` – minimum/maximum valve resistances in the diode law above
 *
 * ### Notes
 * - Set \(R_{\max}\gg R_{\min}\) to suppress backflow while allowing negligible
 *   leakage; \(R_{\min}>0\) adds realistic dissipation during forward flow.
 * - With this formulation, an “ideal” valve would be \(R_{\min}=0\), \(R_{\max}=+\infty\),
 *   but we avoid that for numerical stability, matching the figures.
 */

class PiecewiseCosineChamber : public Block {
 public:
  /**
   * @brief Construct a new BloodVessel object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  PiecewiseCosineChamber(int id, Model *model)
      : Block(id, model, BlockType::piecewise_cosine_chamber, BlockClass::chamber,
              {{"Emax", InputParameter()},
               {"Epass", InputParameter()},
               {"Vrest", InputParameter()},
               {"contract_start", InputParameter()},
               {"relax_start", InputParameter()},
               {"contract_duration", InputParameter()},
               {"relax_duration", InputParameter()}}) {}

  /**
   * @brief Local IDs of the parameters
   *
   */
  enum ParamId {
    EMAX = 0,
    EPASS = 1,
    VREST = 2,
    CSTART = 3,
    RSTART = 4,
    CDUR = 5,
    RDUR = 6
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
   * @brief Number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  TripletsContributions num_triplets{6, 2, 0};

 private:
  double Elas;  // Chamber Elastance

  /**
   * @brief Update the elastance functions which depend on time
   *
   * @param parameters Parameters of the model
   */
  void get_elastance_values(std::vector<double> &parameters);
};

#endif  // SVZERODSOLVER_MODEL_PiecewiseCosineChamber_HPP_
