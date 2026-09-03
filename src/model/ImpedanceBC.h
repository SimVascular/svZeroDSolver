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
 * This boundary condition applies an Olufsen-style discrete convolution using
 * one-cycle flow history:
 *
 * \f[
 * P_{n+1}=P_d+z_0Q_{n+1}+\sum_{m=1}^{N_k-1}z_mQ_{n+1-m}
 * \f]
 *
 * The `z` entries are discrete-time convolution weights. The `m=0` term is
 * treated implicitly through the linear system matrix.
 * History terms (`m>=1`) are evaluated from committed (accepted) timesteps only
 * via a fixed-size ring buffer of one cardiac period.
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
 * `convolution_mode` may be `exact` (default) or `truncated`.
 * In `truncated` mode, only the first `num_kernel_terms` kernel entries are
 * used.
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
   */
  void configure(const std::vector<double>& z, double pd,
                 const std::string& convolution_mode, int num_kernel_terms);

  /**
   * @brief Nonzero triplet counts contributed by this block.
   */
  TripletsContributions num_triplets{2, 0, 0};

 private:
  enum class ConvolutionMode { exact = 0, truncated = 1 };

  void ensure_initialized();
  double lagged_flow(int lag) const;

  std::vector<double> kernel;
  double pd = 0.0;
  ConvolutionMode mode = ConvolutionMode::exact;
  int num_kernel_terms_input = -1;
  bool configured = false;

  // Runtime-discretized state.
  bool initialized = false;
  double dt = -1.0;
  int num_period_steps = -1;
  int num_kernel_terms = -1;
  int head = -1;
  int committed_samples = 0;
  long long num_accepted_steps = 0;
  std::vector<double> flow_history;
};

#endif  // SVZERODSOLVER_MODEL_IMPEDANCEBC_HPP_
