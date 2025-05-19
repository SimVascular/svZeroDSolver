// Test interfacing to svZeroSolver.
// This test mimics an external 3D solver (svSolver/svFSI) interfacing with svZeroDSolver
// The model consists of an RCR BC which acts as a Neumann BC for an external solver
// It mimics two consecutive time steps of an external solver

#include "../LPNSolverInterface/LPNSolverInterface.h" 
#include <iostream>
#include <map>
#include <fstream>
#include <filesystem>
namespace fs = std::filesystem;

//------
// main
//------
//
int main(int argc, char** argv)
{
  LPNSolverInterface interface;

  if (argc != 3) {
    std::runtime_error("Usage: svZeroD_interface_test03 <path_to_svzeroDSolver_build_folder> <path_to_json_file>");
  }
  
  // Load shared library and get interface functions.
  // File extension of the shared library depends on the system
  fs::path build_dir   = argv[1];
  fs::path iface_dir   = build_dir / "src" / "interface";
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
    throw std::runtime_error("Could not find shared libraries " + lib_so.string() + " or " + lib_dylib.string() + " or " + lib_dll.string() + " !");
  }

  // Set up the svZeroD model
  std::string file_name = std::string(argv[2]);
  interface.initialize(file_name);
 
  // Check number of variables and blocks
  if (interface.system_size_ != 3) {
    throw std::runtime_error("interface.system_size_ != 3");
  }
  if (interface.block_names_.size() != 2) {
    throw std::runtime_error("interface.block_names_.size() != 2");
  }

  // Set external time step size (flow solver step size)
  double external_step_size = 0.005;
  interface.set_external_step_size(external_step_size);
  
  // Save the initial condition
  std::vector<double> init_state_y = {-6.2506662304695681e+01, -3.8067539421845140e+04, -3.0504233282976966e+04};
  std::vector<double> init_state_ydot = {-3.0873806830951793e+01, -2.5267653962355386e+05, -2.4894080899699836e+05};

  // Get variable IDs for inlet to RCR block
  std::vector<int> IDs;
  std::string block_name = "RCR"; 
  interface.get_block_node_IDs(block_name, IDs);
  int num_inlet_nodes = IDs[0];
  int num_outlet_nodes = IDs[1+num_inlet_nodes*2];
  if ((num_inlet_nodes != 1) || (num_outlet_nodes != 0)) {
    throw std::runtime_error("Wrong number of inlets/outlets for RCR");
  }
  int rcr_inlet_flow_id = IDs[1];
  int rcr_inlet_pressure_id = IDs[2];
 
  // Update block parameters with current flow from 3D solver
  std::vector<double> new_params(5);
  std::vector<double> params = {-6.2506662041472836e+01, -6.2599344518688739e+01};
  std::vector<double> interface_times = {1.9899999999999796e+00, 1.9949999999999795e+00};
  // Format of new_params for flow/pressure blocks: 
  // [N, time_1, time_2, ..., time_N, value1, value2, ..., value_N]
  // where N is number of time points and value* is flow/pressure
  new_params[0] = 2.0;
  for (int i = 0; i < 2; i++) {
    new_params[1+i] = interface_times[i];
    new_params[3+i] = params[i];
  }
  interface.update_block_params("RCR_coupling", new_params); 
  
  // Set up vectors to run svZeroD simulation
  std::vector<double> solutions(interface.system_size_*interface.num_output_steps_);
  std::vector<double> times(interface.num_output_steps_);
  int error_code = 0;
  
  // Run svZeroD simulation
  interface.update_state(init_state_y, init_state_ydot); 
  interface.run_simulation(interface_times[0], times, solutions, error_code);
  
  // Parse output and calculate mean flow/pressure in aorta and coronary
  int sol_idx = 0;
  double mean_pressure = 0.0;
  for (int tstep = 0; tstep < interface.num_output_steps_; tstep++) {
    for (int state = 0; state < interface.system_size_; state++) {
      sol_idx = interface.system_size_*tstep + state;
      if (state == rcr_inlet_pressure_id) {
        mean_pressure += solutions[sol_idx];
      }
    }
  }
  mean_pressure /= (double)interface.num_output_steps_;
  std::cout <<"Simulation output: " << std::endl;
  std::cout <<"Mean pressure = " << mean_pressure << std::endl;
  
  // Compare mean pressure with pre-computed ("correct") values
  double error_limit = 0.05;
  if (abs(-mean_pressure / 38690.2 - 1.0) > error_limit) {
    throw std::runtime_error("Error in mean pressure at RCR inlet.");
  }

  // Get state vector to prepare for next 3D time step
  interface.return_y(init_state_y);
  interface.return_ydot(init_state_ydot);

  // Update parameters for next 3D time step
  params = {-6.2599344283486374e+01, -6.2630248751732964e+01};
  interface_times = {1.9949999999999795e+00, 1.9999999999999793e+00};
  for (int i = 0; i < 2; i++) {
    new_params[1+i] = interface_times[i];
    new_params[3+i] = params[i];
  }
  interface.update_block_params("RCR_coupling", new_params); 
  
  // Run svZeroD simulation
  interface.update_state(init_state_y, init_state_ydot); 
  interface.run_simulation(interface_times[0], times, solutions, error_code);
  
  // Parse output and calculate mean flow/pressure in aorta and coronary
  sol_idx = 0;
  mean_pressure = 0.0;
  for (int tstep = 0; tstep < interface.num_output_steps_; tstep++) {
    for (int state = 0; state < interface.system_size_; state++) {
      sol_idx = interface.system_size_*tstep + state;
      if (state == rcr_inlet_pressure_id) {
        mean_pressure += solutions[sol_idx];
      }
    }
  }
  mean_pressure /= (double)interface.num_output_steps_;
  std::cout <<"Simulation output: " << std::endl;
  std::cout <<"Mean pressure = " << mean_pressure << std::endl;

  // Compare mean pressure with pre-computed ("correct") values
  if (abs(-mean_pressure / 39911.3 - 1.0) > error_limit) {
    throw std::runtime_error("Error in mean pressure at RCR inlet.");
  }
}
