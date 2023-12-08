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
 * @file svzerodsolver.cpp
 * @brief Main routine of svZeroDSolver
 */
#include <fstream>

#include "Solver.h"

/**
 *
 * @brief svZeroDSolver main routine
 *
 * This is the main routine of the svZeroDSolver. It exectutes the following
 * steps:
 *
 * 1. Read the input file
 * 2. Create the 0D model
 * 3. (Optional) Solve for steady initial condition
 * 4. Run simulation
 * 5. Write output to file
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 * @return Return code
 */
int main(int argc, char* argv[]) {
  DEBUG_MSG("Starting svZeroDSolver");

  // Get input and output file name
  if (argc < 2 || argc > 3) {
    throw std::runtime_error("Usage: svzerodsolver path/to/config.json [optional:path/to/output.csv]");
  }

  std::string input_file_name = argv[1];
  std::string output_file_name;

  if (argc == 3) {
    output_file_name = argv[2];

  } else {
    // If output file is not provided, default is <path to .json>+"output.csv"
    std::size_t end_of_path = input_file_name.rfind("/");

    if (end_of_path == std::string::npos) {
      end_of_path = input_file_name.rfind("\\");  // For Windows paths (?)

      if (end_of_path == std::string::npos) {
        throw std::runtime_error("Error: No output file path provided. Tried to create a default output file but could not find the simulation directory from the input JSON file path.");
      }
    }

    output_file_name = input_file_name.substr(0, end_of_path) + "/output.csv";
    std::cout << "[svzerodsolver] Output will be written to '" << output_file_name << "'." << std::endl;;
  }

  std::ifstream input_file(input_file_name);

  if (!input_file.is_open()) {
    std::cerr << "[svzerodsolver] Error: The input file '" << input_file_name << "' cannot be opened." << std::endl;
    return 1;
  }

  nlohmann::json config;

  try { 
    config = nlohmann::json::parse(input_file);

  } catch (const nlohmann::json::parse_error& e) {
    std::cout << "[svzerodsolver] Error: Parsing the input file '" << input_file_name << "' has failed." << std::endl;
    std::cout << "[svzerodsolver] Details of the parsing error: " << std::endl;
    std::cout << e.what() << std::endl;
    return 1;
  }

  auto solver = Solver(config);
  solver.run();
  solver.write_result_to_csv(output_file_name);

  return 0;
}
