// Test coupling trial/accept semantics for IMPEDANCE BC persistent memory.
//
// The IMPEDANCE block only applies its lagged convolution terms once a full
// cardiac cycle of timesteps has been accepted (see ImpedanceBC::update_time).
// Before that point the block reduces to P = Pd + z[0] * Q, which carries no
// memory at all, so trial/accept semantics are not observable. This test
// therefore warms the block up over one full cycle before asserting anything.
//
// The warm-up also has to impose a *different* flow on each accepted step. The
// convolution is a weighted sum over the one-cycle ring buffer, so a constant
// flow history would produce the same conv_sum before and after a commit and
// the accept step would be unobservable by construction.

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../LPNSolverInterface/LPNSolverInterface.h"
namespace fs = std::filesystem;

namespace {

// Matches svzerod_3Dcoupling_impedance.json: cardiac_period / external step
// size = 0.004 / 0.001, and z.size() == 4.
constexpr double kExternalStepSize = 0.001;
constexpr int kNumPeriodSteps = 4;

// Impose a constant flow `q` over the next external step on the inlet
// coupling block. update_block_params takes a flat FLOW-block payload of
// {num_time_pts, t..., Q...}.
void impose_flow(LPNSolverInterface& interface, double q) {
  std::vector<double> params = {2.0, 0.0, kExternalStepSize, q, q};
  interface.update_block_params("FLOW_COUPLING", params);
}

double max_abs_diff(const std::vector<double>& a,
                    const std::vector<double>& b) {
  double diff = 0.0;
  for (size_t i = 0; i < a.size(); i++) {
    diff = std::max(diff, std::abs(a[i] - b[i]));
  }
  return diff;
}

}  // namespace

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
  interface.set_external_step_size(kExternalStepSize);

  std::vector<double> y(interface.system_size_, 0.0);
  std::vector<double> ydot(interface.system_size_, 0.0);
  interface.return_y(y);
  interface.return_ydot(ydot);

  std::vector<double> t(interface.num_output_steps_, 0.0);
  const size_t soln_size =
      static_cast<size_t>(interface.num_output_steps_) * interface.system_size_;
  std::vector<double> sol_warmup(soln_size, 0.0);
  std::vector<double> sol1(soln_size, 0.0);
  std::vector<double> sol2(soln_size, 0.0);
  std::vector<double> sol3(soln_size, 0.0);

  int error_code = 0;

  // Warm up one full cardiac cycle of accepted steps so the convolution
  // history is live. Each step imposes a distinct flow, and return_y /
  // return_ydot mark the step as accepted (committing the persistent state).
  for (int i = 0; i < kNumPeriodSteps; i++) {
    impose_flow(interface, 1.0 + static_cast<double>(i));
    interface.update_state(y, ydot);
    interface.run_simulation(0.0, t, sol_warmup, error_code);
    if (error_code != 0) {
      throw std::runtime_error("Warm-up run " + std::to_string(i) +
                               " failed.");
    }
    interface.return_y(y);
    interface.return_ydot(ydot);
  }

  // Two trial runs from the same committed state must be identical: the
  // persistent memory is rolled back by update_state, so a rejected trial
  // step must not leak into the next one.
  impose_flow(interface, 9.0);

  interface.update_state(y, ydot);
  interface.run_simulation(0.0, t, sol1, error_code);
  if (error_code != 0) {
    throw std::runtime_error("First trial run failed.");
  }

  interface.update_state(y, ydot);
  interface.run_simulation(0.0, t, sol2, error_code);
  if (error_code != 0) {
    throw std::runtime_error("Second trial run failed.");
  }

  if (max_abs_diff(sol1, sol2) > 1.0e-12) {
    throw std::runtime_error(
        "Repeated trial runs from same committed state are not deterministic "
        "for IMPEDANCE BC.");
  }

  // Accept the trial step, which advances the ring buffer by one sample, then
  // repeat the identical run. Because the history now holds a different set of
  // lagged flows, the result must change.
  std::vector<double> ydot_commit(interface.system_size_, 0.0);
  interface.return_ydot(ydot_commit);

  interface.update_state(y, ydot);
  interface.run_simulation(0.0, t, sol3, error_code);
  if (error_code != 0) {
    throw std::runtime_error("Post-commit run failed.");
  }

  if (max_abs_diff(sol3, sol1) < 1.0e-8) {
    throw std::runtime_error(
        "Persistent state commit had no observable effect for IMPEDANCE BC.");
  }

  std::cout << "test_04 passed" << std::endl;
  return 0;
}
