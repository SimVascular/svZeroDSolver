// Test interfacing to svZeroSolver. 
// This test mimics an external 3D solver (svSolver/svFSI) interfacing with svZeroDSolver
// The model consists of a closed-loop heart model with coronary BCs
// It is run for one time step of the external solver

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
    std::runtime_error("Usage: svZeroD_interface_test01 <path_to_svzeroDSolver_build_folder> <path_to_json_file>");
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
  if (interface.system_size_ != 133) {
    throw std::runtime_error("interface.system_size_ != 133");
  }
  if (interface.block_names_.size() != 50) {
    throw std::runtime_error("interface.block_names_.size() != 50");
  }

  // Set external time step size (flow solver step size)
  double external_step_size = 0.000923;
  interface.set_external_step_size(external_step_size);
  
  // Save the initial condition
  std::vector<double> init_state_y = {220.655, 113.454, 0.146379, 107.558, 0.0840239, 108.31, 0.0917877, 108.966, 0.0539358, 109.893, 0.0997981, 107.152, 0.168397, 109.693, 0.0478851, 108.683, 0.0969178, 106.587, 0.0745793, 111.186, 0.117854, 109.86, 0.063784, 108.403, 0.131471, 110.377, 0.326023, 101.013, 0.127284, 101.488, 0.27798, 110.105, 0.148945, 103.229, 0.14454, 103.893, 0.221119, 104.849, 0.127339, 101.74, 0.156511, 102.527, 0.162979, 103.859, 0.172369, 103.141, 57.563, 1.64141, 54.3487, 1.64141, 0.223534, 1.64141, 0.124233, 1.64141, 0.135591, 1.64141, 0.0763416, 1.64141, 0.151687, 1.64141, 0.253774, 1.64141, 0.0683957, 1.64141, 0.148502, 1.64141, 0.105813, 1.64141, 0.174386, 1.64141, 0.0934222, 1.64141, 0.193053, 1.64141, 0.268681, 1.64141, 0.0993405, 1.64141, 0.211272, 1.64141, 0.11724, 1.64141, 0.112843, 1.64141, 0.175487, 1.64141, 0.0993353, 1.64141, 0.121866, 1.64141, 0.125395, 1.64141, 0.134067, 1.64141, 223.7, 113.546, 223.7, 113.546, 81.4203, -0.00625658, -0.00343448, -0.00367393, -0.00204192, -0.00426901, -0.00667232, -0.00188009, -0.00422262, -0.00273856, -0.00460985, -0.00256085, -0.00507891, -6.67398e-05, 8.96751e-05, 0.0014474, 0.00033538, 0.000384206, 0.000664401, 0.000112356, 0.000156257, 0.000271718, 0.000217284, 35.0055, 1.72547e-12, 45.0271, 68.2839, 555.623, 25.0539, 4.60447, 60.4161, 2.74931e-10, 123.048, 74.8772, 405.595};
  
  std::vector<double> init_state_ydot = {-407.383, 603.025, -0.12541, 586.776, -0.143589, 579.533, -0.143206, 573.381, -0.140919, 563.241, -0.117593, 583.876, -0.149131, 559.217, -0.125706, 567.295, -0.111758, 583.808, -0.152064, 563.677, -0.16802, 563.084, -0.123088, 571.513, -0.201867, 564.857, 1.0169, 633.091, 0.273686, 633.771, 0.040692, 533.643, 0.191667, 575.913, 0.19006, 577.538, 0.223164, 553.187, 0.252782, 625.121, 0.35634, 639.502, 0.327164, 633.898, 0.355656, 632.605, 466.061, -19.3723, 464.294, -19.3723, 0.309113, -19.3723, 0.191591, -19.3723, 0.208691, -19.3723, 0.128617, -19.3723, 0.222605, -19.3723, 0.358526, -19.3723, 0.115032, -19.3723, 0.216634, -19.3723, 0.174236, -19.3723, 0.262547, -19.3723, 0.150376, -19.3723, 0.288608, -19.3723, -0.19198, -19.3723, -0.0722613, -19.3723, -0.0502622, -19.3723, -0.0729683, -19.3723, -0.0668817, -19.3723, -0.0913291, -19.3723, -0.070668, -19.3723, -0.0819415, -19.3723, -0.0760433, -19.3723, -0.085503, -19.3723, -404.441, 515.61, -404.441, 515.61, 662.168, -0.139623, -0.0804254, -0.0861598, -0.0501321, -0.0982909, -0.151062, -0.0461619, -0.0972004, -0.0662369, -0.106878, -0.0616797, -0.116491, -0.0378684, -0.0180084, -0.00754156, -0.0171239, -0.0157096, -0.0189088, -0.0175604, -0.01926, -0.0174746, -0.0195021, 57.5594, -115695, -23.5015, -555.659, -4642.24, 257.474, 24.2288, 555.659, -315843, 345.199, -405.655, -7759.96};

  // Interface blocks flow boundary conditions (neumann boundary conditions for the 3D flow solver)
  std::map<std::string,std::vector<double>> interface_block_params = {
                                            {"inlet_BC_RCR1", {220.655, 220.143}},
                                            {"inlet_BC_lca1", {0.146379, 0.146236}},
                                            {"inlet_BC_lca10", {0.0840239, 0.0838506}},
                                            {"inlet_BC_lca11", {0.0917877, 0.0916064}},
                                            {"inlet_BC_lca12", {0.0539358, 0.0537849}},
                                            {"inlet_BC_lca2", {0.0997981, 0.0996647}},
                                            {"inlet_BC_lca3", {0.168397, 0.168252}},
                                            {"inlet_BC_lca4", {0.0478851, 0.0477658}},
                                            {"inlet_BC_lca5", {0.0969178, 0.0967786}},
                                            {"inlet_BC_lca6", {0.0745793, 0.0743957}},
                                            {"inlet_BC_lca7", {0.117854, 0.117657}},
                                            {"inlet_BC_lca8", {0.063784, 0.0636423}},
                                            {"inlet_BC_lca9", {0.131471, 0.131264}},
                                            {"inlet_BC_rca1", {0.326023, 0.326807}},
                                            {"inlet_BC_rca10", {0.127284, 0.12748}}, 
                                            {"inlet_BC_rca2", {0.27798, 0.278109}},
                                            {"inlet_BC_rca3", {0.148945, 0.149096}},
                                            {"inlet_BC_rca4", {0.14454, 0.144655}},
                                            {"inlet_BC_rca5", {0.221119, 0.221271}},
                                            {"inlet_BC_rca6", {0.127339, 0.127523}},
                                            {"inlet_BC_rca7", {0.156511, 0.15675}},
                                            {"inlet_BC_rca8", {0.162979, 0.163184}}, 
                                            {"inlet_BC_rca9", {0.172369, 0.172628}},
                                            {"outlet_aorta", {223.7, 223.19}}};
  std::vector<double> interface_times = {0.082147, 0.08307};

  // Get variable IDs
  // --- For all interface blocks
  for (const auto block_params : interface_block_params) {
    std::vector<int> IDs;
    std::string block_name = block_params.first;
    double inlet_or_outlet;
    interface.get_block_node_IDs(block_name, IDs);
    // IDs in the above function stores info in the following format:
    // {num inlet nodes, inlet flow[0], inlet pressure[0],..., num outlet nodes, outlet flow[0], outlet pressure[0],...}
    int num_inlet_nodes = IDs[0];
    int num_outlet_nodes = IDs[1+num_inlet_nodes*2];
    if (block_name == "outlet_aorta") {
      if ((num_inlet_nodes != 1) || (num_outlet_nodes != 0)) {
        throw std::runtime_error("Wrong number of inlets/outlets for outlet_aorta");
      }
    } else {
      if ((num_inlet_nodes != 0) || (num_outlet_nodes != 1)) {
        throw std::runtime_error("Wrong number of inlets/outlets for " + block_name);
      }
    }
  }
  // --- For outlet from heart block
  std::vector<int> IDs;
  std::string block_name = "J_heart_outlet"; 
  interface.get_block_node_IDs(block_name, IDs);
  int num_inlet_nodes = IDs[0];
  int num_outlet_nodes = IDs[1+num_inlet_nodes*2];
  if ((num_inlet_nodes != 1) && (num_outlet_nodes != 1)) {
    throw std::runtime_error("Wrong number of inlets/outlets for J_heart_outlet");
  }
  int aortic_inlet_flow_id = IDs[1];
  int aortic_inlet_pressure_id = IDs[2];
  // --- For outlet from coronary
  block_name = "BC_lca1"; 
  interface.get_block_node_IDs(block_name, IDs);
  num_inlet_nodes = IDs[0];
  num_outlet_nodes = IDs[1+num_inlet_nodes*2];
  if ((num_inlet_nodes != 1) && (num_outlet_nodes != 1)) {
    throw std::runtime_error("Wrong number of inlets/outlets for BC_lca1");
  }
  int bc_lca1_outlet_flow_id = IDs[4];
  int bc_lca1_outlet_pressure_id = IDs[5];
 
  // Update block parameters with current flow from 3D solver
  for (const auto block_params : interface_block_params) {
    std::vector<double> new_params(5);
    std::vector<double> params = block_params.second;
    // Format of new_params for flow/pressure blocks: 
    // [N, time_1, time_2, ..., time_N, value1, value2, ..., value_N]
    // where N is number of time points and value* is flow/pressure
    new_params[0] = 2.0;
    for (int i = 0; i < 2; i++) {
      new_params[1+i] = interface_times[i];
      new_params[3+i] = params[i];
    }
    interface.update_block_params(block_params.first, new_params); 
  }
  
  // Set up vectors to run svZeroD simulation
  std::vector<double> solutions(interface.system_size_*interface.num_output_steps_);
  std::vector<double> times(interface.num_output_steps_);
  int error_code = 0;
  
  // Run svZeroD simulation
  interface.update_state(init_state_y, init_state_ydot); 
  interface.run_simulation(0.0, times, solutions, error_code);
  
  // Parse output and calculate mean flow/pressure in aorta and coronary
  int sol_idx = 0;
  double mean_aortic_flow = 0.0;
  double mean_aortic_pressure = 0.0;
  double mean_bc_lca1_outlet_flow = 0.0;
  double mean_bc_lca1_outlet_pressure = 0.0;
  for (int tstep = 0; tstep < interface.num_output_steps_; tstep++) {
    for (int state = 0; state < interface.system_size_; state++) {
      sol_idx = interface.system_size_*tstep + state;
      if (state == aortic_inlet_flow_id) {
        mean_aortic_flow += solutions[sol_idx];
      } else if (state == aortic_inlet_pressure_id) {
        mean_aortic_pressure += solutions[sol_idx];
      } else if (state == bc_lca1_outlet_flow_id) {
        mean_bc_lca1_outlet_flow += solutions[sol_idx];
      } else if (state == bc_lca1_outlet_pressure_id) {
        mean_bc_lca1_outlet_pressure += solutions[sol_idx];
      }
    }
    std::cout << std::endl;
  }
  mean_aortic_flow /= (double)interface.num_output_steps_;
  mean_aortic_pressure /= (double)interface.num_output_steps_;
  mean_bc_lca1_outlet_flow /= (double)interface.num_output_steps_;
  mean_bc_lca1_outlet_pressure /= (double)interface.num_output_steps_;
  
  std::cout <<"Simulation output: " << std::endl;
  std::cout <<"Mean aortic flow = " << mean_aortic_flow << std::endl;
  std::cout <<"Mean aortic pressure = " << mean_aortic_pressure << std::endl;
  std::cout <<"Mean BC_lca1 outlet flow = " << mean_bc_lca1_outlet_flow << std::endl;
  std::cout <<"Mean BC_lca1 outlet pressure = " << mean_bc_lca1_outlet_pressure << std::endl;

  // Compare mean flow/pressure in aorta and coronary with pre-computed ("correct") values
  double error_limit = 0.05;
  std::vector<std::string> wrong_quantities;
  bool is_wrong = false;
  if (abs(mean_aortic_flow / 268.23 - 1.0) > error_limit) {
    is_wrong = true;
    wrong_quantities.push_back("Mean aortic flow");
  }
  if (abs(mean_aortic_pressure / 113.443 - 1.0) > error_limit) {
    is_wrong = true;
    wrong_quantities.push_back("Mean aortic pressure");
  }
  if (abs(mean_bc_lca1_outlet_flow / 0.00755739 - 1.0) > error_limit) {
    is_wrong = true;
    wrong_quantities.push_back("Mean BC_lca1 outlet flow");
  }
  if (abs(mean_bc_lca1_outlet_pressure / 5.88295 - 1.0) > error_limit) {
    is_wrong = true;
    wrong_quantities.push_back("Mean BC_lca1 outlet pressure");
  }
  
  if (is_wrong) {
    std::string error_msg = " ";
    for (int i = 0; i < wrong_quantities.size(); i++) {
      error_msg  = error_msg + wrong_quantities[i] + "; ";
    }
    throw std::runtime_error("Error in the following quantities:" + error_msg);
  }
  
  return 0;
}
