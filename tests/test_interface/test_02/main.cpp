// Test interfacing to svZeroDSolver.
// This test mimics the coupling of svZeroDSolver with an external parameter estimation code (eg. Tulip)
// A coronary model is used and parameters of the BCs are updated and then compared to a reference solution

#include "../LPNSolverInterface/LPNSolverInterface.h" 
#include <iostream>
#include <map>
#include <fstream>
#include <filesystem>
namespace fs = std::filesystem;

//---------------------------------------------------------------------------------------
// Compare mean flow/pressure in aorta and coronary with pre-computed ("correct") values
//---------------------------------------------------------------------------------------
std::string check_simulation_results(double mean_aortic_flow, double mean_aortic_pressure, double mean_bc_lca1_outlet_flow, double mean_bc_lca1_outlet_pressure, double correct_mean_aortic_flow, double correct_mean_aortic_pressure, double correct_mean_bc_lca1_outlet_flow, double correct_mean_bc_lca1_outlet_pressure, bool &is_wrong) {
  double error_limit = 0.01;
  std::vector<std::string> wrong_quantities;
  if (abs(mean_aortic_flow - correct_mean_aortic_flow)/correct_mean_aortic_flow > error_limit) {
    is_wrong = true;
    wrong_quantities.push_back("Mean aortic flow");
  }
  if (abs(mean_aortic_pressure - correct_mean_aortic_pressure)/correct_mean_aortic_pressure > error_limit) {
    is_wrong = true;
    wrong_quantities.push_back("Mean aortic pressure");
  }
  if (abs(mean_bc_lca1_outlet_flow - correct_mean_bc_lca1_outlet_flow)/correct_mean_bc_lca1_outlet_flow > error_limit) {
    is_wrong = true;
    wrong_quantities.push_back("Mean BC_lca1 outlet flow");
  }
  if (abs(mean_bc_lca1_outlet_pressure - correct_mean_bc_lca1_outlet_pressure)/correct_mean_bc_lca1_outlet_pressure > error_limit) {
    is_wrong = true;
    wrong_quantities.push_back("Mean BC_lca1 outlet pressure");
  }
  std::string error_msg = " ";
  if (is_wrong) {
    for (int i = 0; i < wrong_quantities.size(); i++) {
      error_msg = error_msg + wrong_quantities[i] + "; ";
    }
  }
  return error_msg;
}


//------
// main
//------
int main(int argc, char** argv)
{
  LPNSolverInterface interface;
  
  if (argc != 3) {
    std::runtime_error("Usage: svZeroD_interface_test01 <path_to_svZeroDSolver_build_folder> <path_to_json_file>");
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
  if (interface.system_size_ != 253) {
    throw std::runtime_error("interface.system_size_ != 133");
  }
  if (interface.block_names_.size() != 87) {
    throw std::runtime_error("interface.block_names_.size() != 50");
  }

  // Save important variable IDs
  // --- For outlet from heart block
  std::vector<int> IDs;
  std::string block_name = "J_heart_outlet"; 
  interface.get_block_node_IDs(block_name, IDs);
  // IDs in the above function stores info in the following format:
  // {num inlet nodes, inlet flow[0], inlet pressure[0],..., num outlet nodes, outlet flow[0], outlet pressure[0],...}
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
 
  // Save the initial condition
  std::vector<double> init_state_y, init_state_ydot; 
  init_state_y.resize(interface.system_size_);
  init_state_ydot.resize(interface.system_size_);
  interface.return_y(init_state_y);
  interface.return_ydot(init_state_ydot);

  // Set up vectors to run svZeroD simulation
  std::vector<double> solutions(interface.system_size_*interface.num_output_steps_);
  std::vector<double> times(interface.num_output_steps_);
  int error_code = 0;
  
  // Run svZeroD simulation
  interface.update_state(init_state_y, init_state_ydot); 
  interface.run_simulation(0.0, times, solutions, error_code);

  // Parse and print output
  int sol_idx = 0;
  std::vector<double> aortic_flow, aortic_pressure;
  aortic_flow.resize(interface.num_output_steps_);
  aortic_pressure.resize(interface.num_output_steps_);
  double mean_aortic_flow = 0.0;
  double mean_aortic_pressure = 0.0;
  std::vector<double> bc_lca1_outlet_flow, bc_lca1_outlet_pressure;
  bc_lca1_outlet_flow.resize(interface.num_output_steps_);
  bc_lca1_outlet_pressure.resize(interface.num_output_steps_);
  double mean_bc_lca1_outlet_flow = 0.0;
  double mean_bc_lca1_outlet_pressure = 0.0;
  for (int tstep = 0; tstep < interface.num_output_steps_; tstep++) {
    for (int state = 0; state < interface.system_size_; state++) {
      sol_idx = interface.system_size_*tstep + state;
      if (state == aortic_inlet_flow_id) {
        aortic_flow[tstep] = solutions[sol_idx];
        mean_aortic_flow += solutions[sol_idx];
      } else if (state == aortic_inlet_pressure_id) {
        aortic_pressure[tstep] = solutions[sol_idx];
        mean_aortic_pressure += solutions[sol_idx];
      } else if (state == bc_lca1_outlet_flow_id) {
        bc_lca1_outlet_flow[tstep] = solutions[sol_idx];
        mean_bc_lca1_outlet_flow += solutions[sol_idx];
      } else if (state == bc_lca1_outlet_pressure_id) {
        bc_lca1_outlet_pressure[tstep] = solutions[sol_idx];
        mean_bc_lca1_outlet_pressure += solutions[sol_idx];
      }
    }
  }
  mean_aortic_flow /= (double)interface.num_output_steps_;
  mean_aortic_pressure /= (double)interface.num_output_steps_;
  mean_bc_lca1_outlet_flow /= (double)interface.num_output_steps_;
  mean_bc_lca1_outlet_pressure /= (double)interface.num_output_steps_;
  
  std::cout << "First simulation: " << std::endl;
  std::cout << "Mean aortic flow = " << mean_aortic_flow << std::endl;
  std::cout << "Mean aortic pressure = " << mean_aortic_pressure << std::endl;
  std::cout << "Mean BC_lca1 outlet flow = " << mean_bc_lca1_outlet_flow << std::endl;
  std::cout << "Mean BC_lca1 outlet pressure = " << mean_bc_lca1_outlet_pressure << std::endl;

  // Check if outputs are correct
  bool is_wrong = false;
  std::string error_msg;
  error_msg = check_simulation_results(mean_aortic_flow, mean_aortic_pressure, mean_bc_lca1_outlet_flow, mean_bc_lca1_outlet_pressure, 63.3137, 101.139, 0.135942, 3.17525, is_wrong);
  if (is_wrong) {
    throw std::runtime_error("After initial simulation, error in the following quantities: "+error_msg);
  }

  // Save the initial coronary params
  std::map<std::string,std::vector<double>> initial_coronary_params = {
                        {"BC_lca1", {145.272, 236.068, 290.065, 0.00010326, 0.00105039}},
                        {"BC_lca10", {264.502, 429.816, 528.13, 6.513e-05, 0.00066247}},
                        {"BC_lca11", {242.343, 393.807, 483.885, 6.966e-05, 0.00070852}},
                        {"BC_lca12", {434.805, 706.557, 868.172, 4.44e-05, 0.00045194}},
                        {"BC_lca2", {215.568, 350.299, 430.424, 7.617e-05, 0.00077533}},
                        {"BC_lca3", {128.321, 208.522, 256.219, 0.00011358, 0.00115556}},
                        {"BC_lca4", {485.259, 788.546, 968.915, 4.083e-05, 0.00041538}},
                        {"BC_lca5", {220.097, 357.658, 439.467, 7.498e-05, 0.00076305}},
                        {"BC_lca6", {312.826, 508.343, 624.619, 5.727e-05, 0.00058227}},
                        {"BC_lca7", {187.955, 305.427, 375.289, 8.467e-05, 0.00086153}},
                        {"BC_lca8", {353.619, 574.63, 706.069, 5.21e-05, 0.00052984}},
                        {"BC_lca9", {169.536, 275.495, 338.511, 9.166e-05, 0.00093274}},
                        {"BC_rca1", {80.2918, 130.474, 160.318, 0.00017264, 0.00141371}},
                        {"BC_rca10", {218.089, 354.394, 435.457, 8.004e-05, 0.00065544}},
                        {"BC_rca2", {105.209, 170.965, 210.07, 0.00014026, 0.00114837}},
                        {"BC_rca3", {186.143, 302.482, 371.671, 9.038e-05, 0.00074037}},
                        {"BC_rca4", {193.761, 314.861, 386.881, 8.767e-05, 0.00071786}},
                        {"BC_rca5", {124.951, 203.045, 249.489, 0.00012286, 0.0010061}},
                        {"BC_rca6", {218.275, 354.697, 435.829, 7.994e-05, 0.00065505}},
                        {"BC_rca7", {178.053, 289.337, 355.519, 9.357e-05, 0.0007661}},
                        {"BC_rca8", {173.617, 282.127, 346.66, 9.54e-05, 0.00078117}},
                        {"BC_rca9", {162.073, 263.368, 323.61, 0.00010053, 0.00082363}}};

  // Check if correct parameters are read and then update block parameters
  double param_update_factor = 2.0;
  const int num_coronary_params = 5;
  std::vector<double> coronary_params;
  coronary_params.resize(num_coronary_params);
  for (int i = 0; i < interface.block_names_.size(); i++) {
    std::string block_name = interface.block_names_[i];
    if ((block_name.substr(0,6) == "BC_lca") || (block_name.substr(0,6) == "BC_rca")) {
      // Read current block parameters and compare with data above
      interface.read_block_params(block_name, coronary_params);
      auto correct_params = initial_coronary_params[block_name];
      for (int j = 0; j < num_coronary_params; j++) {
        if(abs(coronary_params[j]-correct_params[j])/correct_params[j] > 0.01) {
          throw std::runtime_error("Wrong parameters read from block " + block_name);
        }
      }
      // Update parameters by multiplying Ra, Ram and Rv by param_update_factor
      coronary_params[0] *= param_update_factor; 
      coronary_params[1] *= param_update_factor; 
      coronary_params[2] *= param_update_factor; 
      interface.update_block_params(block_name, coronary_params);
      // Verify that parameters were correctly updated
      interface.read_block_params(block_name, coronary_params);
      for (int j = 0; j < 3; j++) {
        correct_params[j] *= param_update_factor;
        if(abs(coronary_params[j]-correct_params[j])/correct_params[j] > 0.01) {
          throw std::runtime_error("Wrong parameters after update for block " + block_name);
        }
      }
    }
  }
  
  // Run svZeroD simulation with new parameters
  interface.update_state(init_state_y, init_state_ydot); 
  interface.run_simulation(0.0, times, solutions, error_code);

  // Parse and print output
  sol_idx = 0;
  mean_aortic_flow = 0.0;
  mean_aortic_pressure = 0.0;
  mean_bc_lca1_outlet_flow = 0.0;
  mean_bc_lca1_outlet_pressure = 0.0;
  for (int tstep = 0; tstep < interface.num_output_steps_; tstep++) {
    for (int state = 0; state < interface.system_size_; state++) {
      sol_idx = interface.system_size_*tstep + state;
      if (state == aortic_inlet_flow_id) {
        aortic_flow[tstep] = solutions[sol_idx];
        mean_aortic_flow += solutions[sol_idx];
      } else if (state == aortic_inlet_pressure_id) {
        aortic_pressure[tstep] = solutions[sol_idx];
        mean_aortic_pressure += solutions[sol_idx];
      } else if (state == bc_lca1_outlet_flow_id) {
        bc_lca1_outlet_flow[tstep] = solutions[sol_idx];
        mean_bc_lca1_outlet_flow += solutions[sol_idx];
      } else if (state == bc_lca1_outlet_pressure_id) {
        bc_lca1_outlet_pressure[tstep] = solutions[sol_idx];
        mean_bc_lca1_outlet_pressure += solutions[sol_idx];
      }
    }
  }
  mean_aortic_flow /= (double)interface.num_output_steps_;
  mean_aortic_pressure /= (double)interface.num_output_steps_;
  mean_bc_lca1_outlet_flow /= (double)interface.num_output_steps_;
  mean_bc_lca1_outlet_pressure /= (double)interface.num_output_steps_;
  
  std::cout << "After parameter update: " << std::endl;
  std::cout << "Mean aortic flow = " << mean_aortic_flow << std::endl;
  std::cout << "Mean aortic pressure = " << mean_aortic_pressure << std::endl;
  std::cout << "Mean BC_lca1 outlet flow = " << mean_bc_lca1_outlet_flow << std::endl;
  std::cout << "Mean BC_lca1 outlet pressure = " << mean_bc_lca1_outlet_pressure << std::endl;
 
  // Check if outputs are correct
  is_wrong = false;
  error_msg = check_simulation_results(mean_aortic_flow, mean_aortic_pressure, mean_bc_lca1_outlet_flow, mean_bc_lca1_outlet_pressure, 63.2715, 101.232, 0.0691725, 3.17279, is_wrong);
  if (is_wrong) {
    throw std::runtime_error("After parameter update, error in the following quantities: "+error_msg);
  }
  
  // Test restarting simulation by overwriting state with prescribed initial condition and restoring initial parameters
  // Check if correct parameters are read and then update block parameters
  for (int i = 0; i < interface.block_names_.size(); i++) {
    std::string block_name = interface.block_names_[i];
    if ((block_name.substr(0,6) == "BC_lca") || (block_name.substr(0,6) == "BC_rca")) {
      // Read current block parameters and compare with data above
      interface.read_block_params(block_name, coronary_params);
      auto initial_params = initial_coronary_params[block_name];
      // Update parameters by dividing Ra, Ram and Rv by param_update_factor
      coronary_params[0] /= param_update_factor; 
      coronary_params[1] /= param_update_factor; 
      coronary_params[2] /= param_update_factor; 
      interface.update_block_params(block_name, coronary_params);
      // Verify that parameters were correctly updated
      interface.read_block_params(block_name, coronary_params);
      for (int j = 0; j < 3; j++) {
        if(abs(coronary_params[j]-initial_params[j])/initial_params[j] > 0.01) {
          throw std::runtime_error("Wrong parameters after update for block " + block_name);
        }
      }
    }
  }

  // Restore initial state and run simulation
  interface.update_state(init_state_y, init_state_ydot); 
  interface.run_simulation(0.0, times, solutions, error_code);

  // Parse and print output
  sol_idx = 0;
  mean_aortic_flow = 0.0;
  mean_aortic_pressure = 0.0;
  mean_bc_lca1_outlet_flow = 0.0;
  mean_bc_lca1_outlet_pressure = 0.0;
  for (int tstep = 0; tstep < interface.num_output_steps_; tstep++) {
    for (int state = 0; state < interface.system_size_; state++) {
      sol_idx = interface.system_size_*tstep + state;
      if (state == aortic_inlet_flow_id) {
        aortic_flow[tstep] = solutions[sol_idx];
        mean_aortic_flow += solutions[sol_idx];
      } else if (state == aortic_inlet_pressure_id) {
        aortic_pressure[tstep] = solutions[sol_idx];
        mean_aortic_pressure += solutions[sol_idx];
      } else if (state == bc_lca1_outlet_flow_id) {
        bc_lca1_outlet_flow[tstep] = solutions[sol_idx];
        mean_bc_lca1_outlet_flow += solutions[sol_idx];
      } else if (state == bc_lca1_outlet_pressure_id) {
        bc_lca1_outlet_pressure[tstep] = solutions[sol_idx];
        mean_bc_lca1_outlet_pressure += solutions[sol_idx];
      }
    }
  }
  mean_aortic_flow /= (double)interface.num_output_steps_;
  mean_aortic_pressure /= (double)interface.num_output_steps_;
  mean_bc_lca1_outlet_flow /= (double)interface.num_output_steps_;
  mean_bc_lca1_outlet_pressure /= (double)interface.num_output_steps_;
  
  std::cout <<"After restart with same parameters: " << std::endl;
  std::cout <<"Mean aortic flow = " << mean_aortic_flow << std::endl;
  std::cout <<"Mean aortic pressure = " << mean_aortic_pressure << std::endl;
  std::cout <<"Mean BC_lca1 outlet flow = " << mean_bc_lca1_outlet_flow << std::endl;
  std::cout <<"Mean BC_lca1 outlet pressure = " << mean_bc_lca1_outlet_pressure << std::endl;

  // Check if outputs are correct
  is_wrong = false;
  error_msg = check_simulation_results(mean_aortic_flow, mean_aortic_pressure, mean_bc_lca1_outlet_flow, mean_bc_lca1_outlet_pressure, 63.329, 101.158, 0.135971, 3.17321, is_wrong);
  if (is_wrong) {
    throw std::runtime_error("After restart simulation, error in the following quantities: "+error_msg);
  }
  
  return 0;
}
