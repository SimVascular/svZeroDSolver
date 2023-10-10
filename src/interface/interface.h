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

using S = SparseSystem;

/**
 * @brief Interface class for calling svZeroD from external programs
 */

class SolverInterface {
 public:
  SolverInterface(const std::string& input_file_name);
  ~SolverInterface();

  static int problem_id_count_;
  static std::map<int, SolverInterface*> interface_list_;

  int problem_id_ = 0;
  std::string input_file_name_;

  // Parameters for the external solver (the calling program).
  // This is set by the external solver via the interface.
  double external_step_size_ = 0.1;

  // These are read in from the input JSON solver configuration file.
  double time_step_size_ = 0.0;
  int num_time_steps_ = 0;
  double absolute_tolerance_ = 0.0;
  int max_nliter_ = 0;
  int time_step_ = 0.0;
  int save_interval_counter_ = 0;
  int output_interval_ = 0;
  int system_size_ = 0;
  int num_output_steps_ = 0;
  int pts_per_cycle_ = 0;
  bool output_last_cycle_only_ = false;

  std::shared_ptr<Model> model_;
  Integrator integrator_;

  State state_;
  std::vector<double> times_;
  std::vector<State> states_;
};
