

#include "optimize/calibrate.hpp"

// Setting scalar type to double
typedef double T;

int main(int argc, char *argv[]) {
  DEBUG_MSG("Starting svZeroDCalibrator");

  // Get input and output file name
  if (argc != 3) {
    std::cout << "Usage: calibrator path/to/config.json path/to/output.json"
              << std::endl;
    exit(1);
  }

  // Get input and output file names
  std::string input_file = argv[1];
  std::string output_file = argv[2];

  // Parse the configuration
  DEBUG_MSG("Parse configuration");
  std::ifstream ifs(input_file);
  const auto& config = nlohmann::json::parse(ifs);

  auto output_config = OPT::calibrate<T>(config);

  // Write optimized simulation config
  std::ofstream o(output_file);
  o << std::setw(4) << output_config << std::endl;
}
