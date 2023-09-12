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
 * @file solver.hpp
 * @brief SOLVE::Solver source file
 */

#include "algebra/Integrator.h"
#include "algebra/State.h"
#include "helpers/debug.h"
#include "helpers/helpers.h"
#include "io/io.h"
#include "io/csvwriter.h"
#include "io/SimulationParameters.h"
#include "model/Model.h"

#ifndef SVZERODSOLVER_SOLVE_SOLVER_HPP_
#define SVZERODSOLVER_SOLVE_SOLVER_HPP_

namespace solve {

/**
 * @brief Class for running 0D simulations.
 *
 * @tparam T Scalar type (e.g. `float`, `double`)
 */
class Solver {
 public:
  /**
   * @brief Construct a new Solver object
   *
   * @param handler Configuration handler
   */
  Solver(const nlohmann::json& config);

  /**
   * @brief Destroy the Solver object
   *
   */
  ~Solver();

  /**
   * @brief Run the simulation
   *
   */
  void run();

  /**
   * @brief Get the full result as a csv encoded string
   *
   * @return std::string Result
   */
  std::string get_full_result();

  /**
   * @brief Get the result of a single DOF over time
   *
   * @param dof_name Name of the degree-of-freedom
   * @return Eigen::VectorXd Result
   */
  Eigen::VectorXd get_single_result(std::string dof_name);

  /**
   * @brief Get the result of a single DOF averaged over time
   *
   * @param dof_name Name of the degree-of-freedom
   * @return T Result
   */
  double get_single_result_avg(std::string dof_name);

  /**
   * @brief Get the time steps of the result
   *
   * @return std::vector<double>
   */
  std::vector<double> get_times();

  /**
   * @brief Update the parameters of a block
   *
   * @param block_name Name of the block
   * @param new_params New parameters
   */
  void update_block_params(std::string block_name, std::vector<double> new_params);

  /**
   * @brief Write the result to a csv file.
   *
   * @param filename
   */
  void write_result_to_csv(std::string filename);

 private:
  zd_model::Model model;
  io::SimulationParameters simparams;
  std::vector<algebra::State> states;
  std::vector<double> times;
  algebra::State initial_state;

  void sanity_checks();
};


}  

#endif  
