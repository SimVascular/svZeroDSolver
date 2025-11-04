// Test interfacing to svZeroSolver.
// This test mimics an external 3D solver (svSolver/svFSI) interfacing with
// svZeroDSolver The model consists of an RCR BC which acts as a Neumann BC for
// an external solver It mimics two consecutive time steps of an external solver

#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <cstdlib>   // std::getenv

#ifdef _WIN32
  #define NOMINMAX
  #include <windows.h>
#endif

#include "../LPNSolverInterface/LPNSolverInterface.h"
namespace fs = std::filesystem;

static inline void flush_now() { std::cout.flush(); std::cerr.flush(); }
static inline bool debug_on() {
  const char* v = std::getenv("SVZERO_DEBUG");
  return v && *v && std::string(v) != "0";
}

#ifdef _WIN32
static std::string win_err(DWORD e) {
  LPVOID buf = nullptr;
  DWORD n = FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM |
                               FORMAT_MESSAGE_IGNORE_INSERTS,
                           NULL, e, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
                           (LPSTR)&buf, 0, NULL);
  std::string s = (n && buf) ? std::string((LPSTR)buf, n) : "Unknown error";
  if (buf) LocalFree(buf);
  return s;
}

// Make sure Windows can find DLLs in iface_dir without touching global PATH
static void add_dll_search_dir(const fs::path& p) {
  using SetDefaultDllDirectories_t = BOOL (WINAPI*)(DWORD);
  auto setDef = (SetDefaultDllDirectories_t)
      GetProcAddress(GetModuleHandleA("kernel32.dll"), "SetDefaultDllDirectories");
  if (setDef) {
    setDef(LOAD_LIBRARY_SEARCH_DEFAULT_DIRS | LOAD_LIBRARY_SEARCH_USER_DIRS);
    AddDllDirectory(reinterpret_cast<PCWSTR>(std::wstring(p.wstring()).c_str()));
  } else {
    SetDllDirectoryW(std::wstring(p.wstring()).c_str()); // older Windows fallback
  }
}

// LoadLibrary preflight so missing dependencies are surfaced clearly
static void preflight_load(const fs::path& dll) {
  std::wcerr << L"[dbg] Preflight LoadLibraryW: " << dll.wstring() << L"\n";
  HMODULE h = LoadLibraryW(dll.wstring().c_str());
  if (!h) {
    DWORD e = GetLastError();
    std::cerr << "[err] LoadLibraryW failed: " << e << " (" << win_err(e) << ")\n";
    std::cerr << "[hint] Ensure this folder is in DLL search path: "
              << dll.parent_path().string() << "\n";
    flush_now();
    throw std::runtime_error("LoadLibrary preflight failed");
  }
  std::cerr << "[ok ] LoadLibraryW succeeded; FreeLibrary()\n";
  FreeLibrary(h);
  flush_now();
}
#endif

//------
// main
//------
//
int main(int argc, char** argv) {
  // Disable output buffering immediately - critical for Windows CI
  std::setvbuf(stdout, nullptr, _IONBF, 0);
  std::setvbuf(stderr, nullptr, _IONBF, 0);
  std::cout.setf(std::ios::unitbuf);
  std::cerr.setf(std::ios::unitbuf);

  LPNSolverInterface interface;

  if (argc != 3) {
    throw std::runtime_error(
        "Usage: svZeroD_interface_test03 <path_to_svzeroDSolver_build_folder> "
        "<path_to_json_file>");
  }

  // Load shared library and get interface functions.
  // File extension of the shared library depends on the system
  fs::path build_dir = argv[1];
  fs::path iface_dir = build_dir / "src" / "interface";
  fs::path lib_so = iface_dir / "libsvzero_interface.so";
  fs::path lib_dylib = iface_dir / "libsvzero_interface.dylib";
#ifdef _WIN32
  fs::path lib_dll1 = iface_dir / "libsvzero_interface.dll"; // with prefix (our CMake patch)
  fs::path lib_dll2 = iface_dir / "svzero_interface.dll";    // fallback without prefix
#endif

  std::string lib_to_load;
#ifdef _WIN32
  if (fs::exists(lib_dll1)) lib_to_load = lib_dll1.string();
  else if (fs::exists(lib_dll2)) lib_to_load = lib_dll2.string();
#else
  if (fs::exists(lib_so)) lib_to_load = lib_so.string();
  else if (fs::exists(lib_dylib)) lib_to_load = lib_dylib.string();
#endif
  if (lib_to_load.empty()) {
    std::cerr << "[err] Could not find shared library in " << iface_dir.string() << "\n";
    flush_now();
    throw std::runtime_error("Could not find shared library!");
  }

#ifdef _WIN32
  add_dll_search_dir(iface_dir);
  preflight_load(fs::path(lib_to_load));
#endif

  interface.load_library(lib_to_load);

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
  std::vector<double> init_state_y = {-6.2506662304695681e+01,
                                      -3.8067539421845140e+04,
                                      -3.0504233282976966e+04};
  std::vector<double> init_state_ydot = {-3.0873806830951793e+01,
                                         -2.5267653962355386e+05,
                                         -2.4894080899699836e+05};

  // Get variable IDs for inlet to RCR block
  std::vector<int> IDs;
  std::string block_name = "RCR";
  interface.get_block_node_IDs(block_name, IDs);
  int num_inlet_nodes = IDs[0];
  int num_outlet_nodes = IDs[1 + num_inlet_nodes * 2];
  if ((num_inlet_nodes != 1) || (num_outlet_nodes != 0)) {
    throw std::runtime_error("Wrong number of inlets/outlets for RCR");
  }
  int rcr_inlet_flow_id = IDs[1];
  int rcr_inlet_pressure_id = IDs[2];

  // Update block parameters with current flow from 3D solver
  std::vector<double> new_params(5);
  std::vector<double> params = {-6.2506662041472836e+01,
                                -6.2599344518688739e+01};
  std::vector<double> interface_times = {1.9899999999999796e+00,
                                         1.9949999999999795e+00};
  // Format of new_params for flow/pressure blocks:
  // [N, time_1, time_2, ..., time_N, value1, value2, ..., value_N]
  // where N is number of time points and value* is flow/pressure
  new_params[0] = 2.0;
  for (int i = 0; i < 2; i++) {
    new_params[1 + i] = interface_times[i];
    new_params[3 + i] = params[i];
  }
  interface.update_block_params("RCR_coupling", new_params);

  // Set up vectors to run svZeroD simulation
  std::vector<double> solutions(interface.system_size_ *
                                interface.num_output_steps_);
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
      sol_idx = interface.system_size_ * tstep + state;
      if (state == rcr_inlet_pressure_id) {
        mean_pressure += solutions[sol_idx];
      }
    }
  }
  mean_pressure /= (double)interface.num_output_steps_;
  std::cout << "Simulation output: " << std::endl;
  std::cout << "Mean pressure = " << mean_pressure << std::endl;

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
    new_params[1 + i] = interface_times[i];
    new_params[3 + i] = params[i];
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
      sol_idx = interface.system_size_ * tstep + state;
      if (state == rcr_inlet_pressure_id) {
        mean_pressure += solutions[sol_idx];
      }
    }
  }
  mean_pressure /= (double)interface.num_output_steps_;
  std::cout << "Simulation output: " << std::endl;
  std::cout << "Mean pressure = " << mean_pressure << std::endl;

  // Compare mean pressure with pre-computed ("correct") values
  if (abs(-mean_pressure / 39911.3 - 1.0) > error_limit) {
    throw std::runtime_error("Error in mean pressure at RCR inlet.");
  }
}
