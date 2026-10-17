// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file ImpedanceBC.h
 * @brief model::ImpedanceBC source file
 */
#ifndef SVZERODSOLVER_MODEL_IMPEDANCEBC_HPP_
#define SVZERODSOLVER_MODEL_IMPEDANCEBC_HPP_

#include <string>
#include <vector>

#include "Block.h"
#include "SparseSystem.h"

/**
 * @brief Structured-tree impedance boundary condition with periodic memory.
 *
 * Imposes the outlet pressure as a discrete convolution of the impedance
 * kernel with one cardiac cycle of flow history, following the structured-tree
 * outflow condition of \cite olufsen2000.
 *
 * ### Governing equations
 *
 * \f[
 * P_{n+1}=P_d+z_0Q_{n+1}+\sum_{m=1}^{N_k-1}z_mQ_{n+1-m}
 * \f]
 *
 * The \f$z\f$ entries are discrete-time convolution weights. The \f$m=0\f$
 * term depends on the current flow, so it enters the system matrix directly;
 * the history terms (\f$m\geq 1\f$) are evaluated from accepted timesteps only,
 * held in a fixed-size ring buffer spanning one cardiac cycle.
 *
 * The ring buffer is zero-initialized, so the history terms are omitted
 * entirely until one full cycle has been accepted, after which they switch on
 * at full magnitude. This reproduces the reference implementation, which gates
 * the same sum on `t > Period`, and makes the pressure discontinuous at the
 * start of the second cycle. The first cycle is therefore not a meaningful
 * solution and should be discarded.
 *
 * ### Local contributions
 *
 * \f[
 * \mathbf{y}^{e}=\left[\begin{array}{ll}P^{e} & Q^{e}\end{array}\right]^{T}
 * \f]
 *
 * \f[
 * \mathbf{F}^{e}=\left[\begin{array}{cc}
 * 1 & -z_{0}
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{c}^{e}=\left[\begin{array}{c}
 * -\left(P_{d}+\sum_{m=1}^{N_k-1}z_mQ_{n+1-m}\right)
 * \end{array}\right]
 * \f]
 *
 * ### Parameters
 *
 * * `z`: Time-domain impedance kernel. Its length must equal one cardiac cycle
 *   in time steps (see the timing requirements below).
 * * `Pd`: Distal/reference pressure. Optional, defaults to `0.0`.
 * * `convolution_mode`: `exact` (default) or `truncated`. In `truncated` mode
 *   only the first `num_kernel_terms` kernel entries are used, which costs
 *   \f$O(N_k)\f$ per accepted timestep instead of one full cycle.
 * * `num_kernel_terms`: Number of kernel terms. Required when
 *   `convolution_mode` is `truncated`, ignored otherwise.
 *
 * ### Usage in json configuration file
 *
 *     "boundary_conditions": [
 *         {
 *             "bc_name": "OUT",
 *             "bc_type": "IMPEDANCE",
 *             "bc_values": {
 *                 "Pd": 0.0,
 *                 "z": [ ... ],
 *                 "convolution_mode": "exact",
 *                 "num_kernel_terms": 128
 *             }
 *         }
 *     ]
 *
 * ### Simulation timing requirements
 *
 * The kernel is discretized at the solver time step, which constrains
 * `simulation_parameters`:
 *
 * * `cardiac_period / dt` must equal `z.size()` in all cases.
 * * Non-coupled configurations must set
 *   `number_of_time_pts_per_cardiac_cycle = z.size() + 1`.
 * * 3D-coupled configurations must keep `number_of_time_pts = 2` for the
 *   external interface step. If `number_of_time_pts_per_cardiac_cycle` is
 *   omitted, the cycle discretization is inferred from `cardiac_period / dt`.
 * * `steady_initial` must be set to `false`. The steady initial condition is
 *   solved with a time step of `cardiac_period / 10`, which does not match the
 *   kernel discretization and has no meaningful interpretation for this
 *   boundary condition. Since `steady_initial` defaults to `true`, it must be
 *   set explicitly; otherwise the solver reports an error and exits.
 */
class ImpedanceBC : public Block {
 public:
  /**
   * @brief Construct a new ImpedanceBC object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  ImpedanceBC(int id, Model* model)
      : Block(id, model, BlockType::impedance_bc, BlockClass::boundary_condition,
              {{"z", InputParameter(false, false, false)},
               {"Pd", InputParameter(true, false, false)},
               {"convolution_mode", InputParameter(true, false, false)},
               {"num_kernel_terms", InputParameter(true, false, false)}}) {}

  /**
   * @brief Set up the algebraic degree of freedom for the impedance BC.
   *
   * @param dofhandler Global degree-of-freedom handler
   */
  void setup_dofs(DOFHandler& dofhandler) override;

  /**
   * @brief Assemble constant Jacobian contributions.
   *
   * @param system Global sparse system
   * @param parameters Current parameter values
   */
  void update_constant(SparseSystem& system,
                       std::vector<double>& parameters) override;

  /**
   * @brief Assemble time-dependent residual contributions.
   *
   * @param system Global sparse system
   * @param parameters Current parameter values
   */
  void update_time(SparseSystem& system,
                   std::vector<double>& parameters) override;

  /**
   * @brief Commit the accepted outlet flow for the next convolution update.
   *
   * @param y Accepted global state vector
   */
  void accept_timestep(
      const Eigen::Matrix<double, Eigen::Dynamic, 1>& y) override;

  /**
   * @brief Serialize the convolution history needed for restart/coupling.
   *
   * @return Persistent block state
   */
  std::vector<double> get_persistent_state() const override;

  /**
   * @brief Restore a previously serialized convolution history state.
   *
   * @param state Persistent block state
   */
  void set_persistent_state(const std::vector<double>& state) override;

  /**
   * @brief Configure the impedance kernel and convolution controls.
   *
   * @param z Time-domain impedance kernel
   * @param pd Distal/reference pressure
   * @param convolution_mode Convolution mode (`exact` or `truncated`)
   * @param num_kernel_terms Number of kernel terms for truncated mode
   * @param points_per_cycle Sample points spanning one cardiac cycle
   * (`simulation_parameters.number_of_time_pts_per_cardiac_cycle`). Pass 0 to
   * infer the cycle discretization from `cardiac_period / dt` instead.
   */
  void configure(const std::vector<double>& z, double pd,
                 const std::string& convolution_mode, int num_kernel_terms,
                 int points_per_cycle);

  /**
   * @brief Nonzero triplet counts contributed by this block.
   */
  TripletsContributions num_triplets{2, 0, 0};

 private:
  enum class ConvolutionMode { exact = 0, truncated = 1 };

  /// @brief Scalar fields preceding the flow history in the persistent state.
  /// Layout: [dt, num_period_steps, num_kernel_terms, head,
  ///          num_accepted_steps, flow_history...]
  static constexpr size_t kNumStateHeaderFields = 5;

  void ensure_initialized();
  double lagged_flow(int lag) const;

  std::vector<double> kernel;
  double pd = 0.0;
  ConvolutionMode mode = ConvolutionMode::exact;
  int num_kernel_terms_input = -1;
  /// @brief One-cycle sample count from the configuration; 0 means infer the
  /// cycle discretization from `cardiac_period / dt`.
  int points_per_cycle_input = 0;
  bool configured = false;

  // Runtime-discretized state.
  bool initialized = false;
  double dt = -1.0;
  int num_period_steps = -1;
  int num_kernel_terms = -1;
  int head = -1;
  /// @brief Number of timesteps accepted since (re)initialization. Also tells
  /// how much of the ring buffer holds real data: min(this, num_period_steps).
  long long num_accepted_steps = 0;
  std::vector<double> flow_history;
};

#endif  // SVZERODSOLVER_MODEL_IMPEDANCEBC_HPP_
