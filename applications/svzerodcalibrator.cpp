// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause
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
