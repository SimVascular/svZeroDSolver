// Test coupling trial/accept semantics for IMPEDANCE BC persistent memory.

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../LPNSolverInterface/LPNSolverInterface.h"
namespace fs = std::filesystem;

int main(int argc, char** argv) {
  LPNSolverInterface interface;

  if (argc != 3) {
    throw std::runtime_error(
        "Usage: svZeroD_interface_test04 <path_to_svzeroDSolver_build_folder> "
        "<path_to_json_file>");
  }

  fs::path build_dir = argv[1];
  fs::path iface_dir = build_dir / "src" / "interface";
  fs::path lib_so = iface_dir / "libsvzero_interface.so";
  fs::path lib_dylib = iface_dir / "libsvzero_interface.dylib";
  fs::path lib_dll = iface_dir / "libsvzero_interface.dll";

  if (fs::exists(lib_so)) {
    interface.load_library(lib_so.string());
  } else if (fs::exists(lib_dylib)) {
    interface.load_library(lib_dylib.string());
  } else if (fs::exists(lib_dll)) {
    interface.load_library(lib_dll.string());
  } else {
    throw std::runtime_error("Could not find shared libraries " +
                             lib_so.string() + " or " + lib_dylib.string() +
                             " or " + lib_dll.string() + " !");
  }

  interface.initialize(std::string(argv[2]));
  interface.set_external_step_size(0.001);

  std::vector<double> y0(interface.system_size_, 0.0);
  std::vector<double> ydot0(interface.system_size_, 0.0);
  interface.return_y(y0);
  interface.return_ydot(ydot0);

  std::vector<double> params = {2.0, 0.0, 0.001, 1.0, 2.0};
  interface.update_block_params("FLOW_COUPLING", params);

  std::vector<double> t(interface.num_output_steps_, 0.0);
  std::vector<double> sol1(interface.num_output_steps_ * interface.system_size_,
                           0.0);
  std::vector<double> sol2(interface.num_output_steps_ * interface.system_size_,
                           0.0);
  std::vector<double> sol3(interface.num_output_steps_ * interface.system_size_,
                           0.0);

  int error_code = 0;

  interface.update_state(y0, ydot0);
  interface.run_simulation(0.0, t, sol1, error_code);
  if (error_code != 0) {
    throw std::runtime_error("First trial run failed.");
  }

  interface.update_state(y0, ydot0);
  interface.run_simulation(0.0, t, sol2, error_code);
  if (error_code != 0) {
    throw std::runtime_error("Second trial run failed.");
  }

  double max_trial_diff = 0.0;
  for (size_t i = 0; i < sol1.size(); i++) {
    max_trial_diff = std::max(max_trial_diff, std::abs(sol1[i] - sol2[i]));
  }
  if (max_trial_diff > 1.0e-12) {
    throw std::runtime_error(
        "Repeated trial runs from same committed state are not deterministic "
        "for IMPEDANCE BC.");
  }

  // Mark accepted state (and persistent memory) using current interface state.
  std::vector<double> ydot_commit(interface.system_size_, 0.0);
  interface.return_ydot(ydot_commit);

  interface.update_state(y0, ydot0);
  interface.run_simulation(0.0, t, sol3, error_code);
  if (error_code != 0) {
    throw std::runtime_error("Post-commit run failed.");
  }

  double max_commit_diff = 0.0;
  for (size_t i = 0; i < sol1.size(); i++) {
    max_commit_diff = std::max(max_commit_diff, std::abs(sol3[i] - sol1[i]));
  }
  if (max_commit_diff < 1.0e-8) {
    throw std::runtime_error(
        "Persistent state commit had no observable effect for IMPEDANCE BC.");
  }

  std::cout << "test_04 passed" << std::endl;
  return 0;
}
