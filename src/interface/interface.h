// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file interface.h
 * @brief svZeroDSolver callable interface.
 */

#include <map>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

#include "Integrator.h"
#include "Model.h"
#include "SparseSystem.h"
#include "State.h"
#include "csv_writer.h"
#include "debug.h"

/**
 * @brief Interface class for calling svZeroD from external programs
 */
class SolverInterface {
 public:
  /**
   * @brief Construct a new interface object
   * @param input_file_name The 0D JSON file which specifies the model
   */
  SolverInterface(const std::string& input_file_name);

  /**
   * @brief Destroy the interface object
   */
  ~SolverInterface();

  /**
   * @brief Counter for the number of interfaces
   */
  static int problem_id_count_;
  /**
   * @brief List of interfaces
   */
  static std::map<int, SolverInterface*> interface_list_;

  /**
   * @brief ID of current interface
   */
  int problem_id_ = 0;

  /**
   * @brief 0D input (JSON) file
   */
  std::string input_file_name_;

  /**
   * @brief Time step size of the external program
   *
   * This is required for coupling with a 3D solver
   */
  double external_step_size_ = 0.1;

  // These are read in from the input JSON solver configuration file.
  /**
   * @brief 0D time step size
   */
  double time_step_size_ = 0.0;

  /**
   * @brief Spectral radius of generalized alpha integrator
   */
  double rho_infty_ = 0.0;

  /**
   * @brief Number of 0D time steps
   */
  int num_time_steps_ = 0;
  /**
   * @brief Convergence tolerance for the 0D model
   */
  double absolute_tolerance_ = 0.0;
  /**
   * @brief Maximum number of non-linear iterations
   */
  int max_nliter_ = 0;
  /**
   * @brief Current time step
   */
  int time_step_ = 0.0;
  /**
   * @brief The size of the 0D system
   */
  int system_size_ = 0;
  /**
   * @brief The number of steps to output
   */
  int num_output_steps_ = 0;
  /**
   * @brief Number of time steps per cycle
   */
  int pts_per_cycle_ = 0;
  /**
   * @brief Output results from last cycle only?
   */
  bool output_last_cycle_only_ = false;

  /**
   * @brief The current 0D model object
   */
  std::shared_ptr<Model> model_;
  /**
   * @brief The current 0D integrator object
   */
  Integrator integrator_;

  /**
   * @brief The current 0D state vector
   */
  State state_;
  /**
   * @brief Vector to store solution times
   */
  std::vector<double> times_;
  /**
   * @brief Vector to store solution states
   */
  std::vector<State> states_;
};
