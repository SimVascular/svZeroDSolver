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
 * @file svzerodcalibrator.cpp
 * @brief Main routine for svZeroDCalibrator
 */
#include "calibrate.h"

int main(int argc, char* argv[]) {
  DEBUG_MSG("Starting svZeroDCalibrator");

  // Get input and output file name
  if (argc != 3) {
    std::cout
        << "Usage: svzerodcalibrator path/to/config.json path/to/output.json"
        << std::endl;
    return 1;
  }

  // Get input and output file names
  std::string input_file_name = argv[1];
  std::string output_file_name = argv[2];

  // Parse the configuration
  DEBUG_MSG("Parse configuration");
  std::ifstream input_file(input_file_name);

  if (!input_file.is_open()) {
    std::cerr << "[svzerodcalibrator] Error: The input file '" << input_file_name << "' cannot be opened." << std::endl;
    return 1;
  }

  const auto& config = nlohmann::json::parse(input_file);
  nlohmann::json output_config;

  try {
    output_config = calibrate(config);
  } catch (const nlohmann::json::parse_error& e) {
    std::cerr << "[svzerodcalibrator] Error: The input file '" << input_file_name
    << "' does not have the parameters needed by the calibrate program." << std::endl;
    return 1;
  }

  // Write optimized simulation config
  std::ofstream out_file(output_file_name);

  if (!out_file.is_open()) {
    std::cerr << "[svzerodcalibrator] Error: The output file '" << output_file_name << "' cannot be opened." << std::endl;
    return 1;
  }

  out_file << std::setw(4) << output_config << std::endl;
}
