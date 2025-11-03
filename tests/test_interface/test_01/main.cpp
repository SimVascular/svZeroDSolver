// Test interfacing to svZeroSolver.
// This test mimics an external 3D solver (svSolver/svFSI) interfacing with
// svZeroDSolver. The model consists of a closed-loop heart model with coronary
// BCs. It is run for one time step of the external solver.

#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <cstdlib>   // std::getenv
#include <thread>
#include <chrono>
#include <atomic>

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

// Simple heartbeat printer so we see progress while inside long calls
struct Heartbeat {
  std::atomic<bool> stop{false};
  std::thread th;
  Heartbeat() {
    if (debug_on()) {
      th = std::thread([this]{
        using namespace std::chrono_literals;
        while (!stop.load()) {
          std::cerr << ".";
          flush_now();
          std::this_thread::sleep_for(1s);
        }
      });
    }
  }
  ~Heartbeat() {
    stop.store(true);
    if (th.joinable()) th.join();
  }
};

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

int main(int argc, char** argv) {
  try {
    if (argc != 3) {
      throw std::runtime_error(
        "Usage: svZeroD_interface_test01 <path_to_svzeroDSolver_build_folder> <path_to_json_file>");
    }

    std::cout << "[dbg] CWD            : " << fs::current_path().string() << "\n";
    std::cout << "[dbg] argv[1] build  : " << argv[1] << "\n";
    std::cout << "[dbg] argv[2] json   : " << argv[2] << "\n";
#ifdef _WIN32
    if (debug_on()) {
      std::cout << "[dbg] (win) initial PATH (first 4):\n";
      int shown = 0;
      for (const auto& part : std::string(std::getenv("PATH") ? std::getenv("PATH") : "").substr(0, 8192)) {
        // just ensure we flush something; detailed PATH dump can be too long
        (void)part; break;
      }
      std::cout << "      (suppressed full PATH dump)\n";
    }
#endif
    flush_now();

    fs::path build_dir = fs::path(argv[1]);
    fs::path iface_dir = build_dir / "src" / "interface";
    fs::path lib_so    = iface_dir / "libsvzero_interface.so";
    fs::path lib_dylib = iface_dir / "libsvzero_interface.dylib";
#ifdef _WIN32
    fs::path lib_dll1  = iface_dir / "libsvzero_interface.dll"; // with prefix (our CMake patch)
    fs::path lib_dll2  = iface_dir / "svzero_interface.dll";    // fallback without prefix
#endif

    std::cout << "[dbg] iface_dir       : " << iface_dir.string() << "\n";
    std::cout << "[dbg] exists(iface)   : " << (fs::exists(iface_dir) ? "yes" : "no") << "\n";
#ifndef _WIN32
    std::cout << "[dbg] so/dylib cand.  : " << lib_so.string() << " | " << lib_dylib.string() << "\n";
#else
    std::cout << "[dbg] dll candidates  : " << lib_dll1.string() << " | " << lib_dll2.string() << "\n";
#endif
    flush_now();

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
      return 2;
    }
    std::cout << "[dbg] using library   : " << lib_to_load << "\n";
    flush_now();

#ifdef _WIN32
    add_dll_search_dir(iface_dir);
    preflight_load(fs::path(lib_to_load));
#endif

    LPNSolverInterface interface;

    std::cout << "[step] load_library...\n"; flush_now();
    interface.load_library(lib_to_load);
    std::cout << "[ok  ] load_library\n"; flush_now();

    // Set up the svZeroD model
    const std::string file_name = std::string(argv[2]);
    std::cout << "[step] initialize: " << file_name << "\n"; flush_now();
    if (!fs::exists(file_name)) {
      std::cerr << "[err] JSON file does not exist: " << file_name << "\n";
      flush_now();
      return 3;
    }
    interface.initialize(file_name);
    std::cout << "[ok  ] initialize\n"; flush_now();

    // Check number of variables and blocks
    std::cout << "[dbg] system_size_    : " << interface.system_size_ << "\n";
    std::cout << "[dbg] block_names_.sz : " << interface.block_names_.size() << "\n";
    flush_now();
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
    std::vector<double> init_state_y = {
        220.655,     113.454,      0.146379,    107.558,     0.0840239,
        108.31,      0.0917877,    108.966,     0.0539358,   109.893,
        0.0997981,   107.152,      0.168397,    109.693,     0.0478851,
        108.683,     0.0969178,    106.587,     0.0745793,   111.186,
        0.117854,    109.86,       0.063784,    108.403,     0.131471,
        110.377,     0.326023,     101.013,     0.127284,    101.488,
        0.27798,     110.105,      0.148945,    103.229,     0.14454,
        103.893,     0.221119,     104.849,     0.127339,    101.74,
        0.156511,    102.527,      0.162979,    103.859,     0.172369,
        103.141,     57.563,       1.64141,     54.3487,     1.64141,
        0.223534,    1.64141,      0.124233,    1.64141,     0.135591,
        1.64141,     0.0763416,    1.64141,     0.151687,    1.64141,
        0.253774,    1.64141,      0.0683957,   1.64141,     0.148502,
        1.64141,     0.105813,     1.64141,     0.174386,    1.64141,
        0.0934222,   1.64141,      0.193053,    1.64141,     0.268681,
        1.64141,     0.0993405,    1.64141,     0.211272,    1.64141,
        0.11724,     1.64141,      0.112843,    1.64141,     0.175487,
        1.64141,     0.0993353,    1.64141,     0.121866,    1.64141,
        0.125395,    1.64141,      0.134067,    1.64141,     223.7,
        113.546,     223.7,        113.546,     81.4203,     -0.00625658,
        -0.00343448, -0.00367393,  -0.00204192, -0.00426901, -0.00667232,
        -0.00188009, -0.00422262,  -0.00273856, -0.00460985, -0.00256085,
        -0.00507891, -6.67398e-05, 8.96751e-05, 0.0014474,   0.00033538,
        0.000384206, 0.000664401,  0.000112356, 0.000156257, 0.000271718,
        0.000217284, 35.0055,      1.72547e-12, 45.0271,     68.2839,
        555.623,     25.0539,      4.60447,     60.4161,     2.74931e-10,
        123.048,     74.8772,      405.595};

    std::vector<double> init_state_ydot = {
        -407.383,   603.025,    -0.12541,   586.776,    -0.143589,  579.533,
        -0.143206,  573.381,    -0.140919,  563.241,    -0.117593,  583.876,
        -0.149131,  559.217,    -0.125706,  567.295,    -0.111758,  583.808,
        -0.152064,  563.677,    -0.16802,   563.084,    -0.123088,  571.513,
        -0.201867,  564.857,    1.0169,     633.091,    0.273686,   633.771,
        0.040692,   533.643,    0.191667,   575.913,    0.19006,    577.538,
        0.223164,   553.187,    0.252782,   625.121,    0.35634,    639.502,
        0.327164,   633.898,    0.355656,   632.605,    466.061,    -19.3723,
        464.294,    -19.3723,   0.309113,   -19.3723,   0.191591,   -19.3723,
        0.208691,   -19.3723,   0.128617,   -19.3723,   0.222605,   -19.3723,
        0.358526,   -19.3723,   0.115032,   -19.3723,   0.216634,   -19.3723,
        0.174236,   -19.3723,   0.262547,   -19.3723,   0.150376,   -19.3723,
        0.288608,   -19.3723,   -0.19198,   -19.3723,   -0.0722613, -19.3723,
        -0.0502622, -19.3723,   -0.0729683, -19.3723,   -0.0668817, -19.3723,
        -0.0913291, -19.3723,   -0.070668,  -19.3723,   -0.0819415, -19.3723,
        -0.0760433, -19.3723,   -0.085503,  -19.3723,   -404.441,   515.61,
        -404.441,   515.61,     662.168,    -0.139623,  -0.0804254, -0.0861598,
        -0.0501321, -0.0982909, -0.151062,  -0.0461619, -0.0972004, -0.0662369,
        -0.106878,  -0.0616797, -0.116491,  -0.0378684, -0.0180084, -0.00754156,
        -0.0171239, -0.0157096, -0.0189088, -0.0175604, -0.01926,   -0.0174746,
        -0.0195021, 57.5594,    -115695,    -23.5015,   -555.659,   -4642.24,
        257.474,    24.2288,    555.659,    -315843,    345.199,    -405.655,
        -7759.96};

    // Interface blocks flow boundary conditions (neumann boundary conditions for the 3D flow solver)
    std::map<std::string, std::vector<double>> interface_block_params = {
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

    std::cout << "[step] query block node IDs...\n"; flush_now();
    // Get variable IDs for all interface blocks
    for (const auto& block_params : interface_block_params) {
      std::vector<int> IDs;
      const std::string& block_name = block_params.first;
      std::cout << "[dbg] Processing block: " << block_name << "\n"; flush_now();
      interface.get_block_node_IDs(block_name, IDs);
      if (IDs.size() < 2) throw std::runtime_error("IDs too small for " + block_name);
      int num_inlet_nodes  = IDs[0];
      int num_outlet_nodes = IDs[1 + num_inlet_nodes * 2];
      std::cout << "[dbg] " << block_name << ": inlets=" << num_inlet_nodes 
          << ", outlets=" << num_outlet_nodes << "\n"; flush_now();
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
    // Outlet from heart block
    std::vector<int> IDs;
    std::string block_name = "J_heart_outlet";
    std::cout << "[dbg] Processing heart outlet block\n"; flush_now();
    interface.get_block_node_IDs(block_name, IDs);
    int num_inlet_nodes = IDs[0];
    int num_outlet_nodes = IDs[1 + num_inlet_nodes * 2];
    std::cout << "[dbg] Heart outlet: inlets=" << num_inlet_nodes 
          << ", outlets=" << num_outlet_nodes << "\n"; flush_now();
    if ((num_inlet_nodes != 1) && (num_outlet_nodes != 1)) {
      throw std::runtime_error("Wrong number of inlets/outlets for J_heart_outlet");
    }
    int aortic_inlet_flow_id = IDs[1];
    int aortic_inlet_pressure_id = IDs[2];
    std::cout << "[dbg] Aortic IDs: flow=" << aortic_inlet_flow_id 
          << ", pressure=" << aortic_inlet_pressure_id << "\n"; flush_now();
    // Outlet from coronary
    block_name = "BC_lca1";
    std::cout << "[dbg] Processing coronary block\n"; flush_now();
    interface.get_block_node_IDs(block_name, IDs);
    num_inlet_nodes = IDs[0];
    num_outlet_nodes = IDs[1 + num_inlet_nodes * 2];
    std::cout << "[dbg] Coronary: inlets=" << num_inlet_nodes 
          << ", outlets=" << num_outlet_nodes << "\n"; flush_now();
    if ((num_inlet_nodes != 1) && (num_outlet_nodes != 1)) {
      throw std::runtime_error("Wrong number of inlets/outlets for BC_lca1");
    }
    int bc_lca1_outlet_flow_id = IDs[4];
    int bc_lca1_outlet_pressure_id = IDs[5];
    std::cout << "[dbg] Coronary IDs: flow=" << bc_lca1_outlet_flow_id 
          << ", pressure=" << bc_lca1_outlet_pressure_id << "\n"; flush_now();
    std::cout << "[ok  ] block node IDs\n"; flush_now();

    // Update block parameters with current flow from 3D solver
    std::cout << "[step] update_block_params...\n"; flush_now();
    for (const auto& block_params : interface_block_params) {
      std::vector<double> new_params(5);
      const std::vector<double>& params = block_params.second;
      new_params[0] = 2.0;
      for (int i = 0; i < 2; i++) {
        new_params[1 + i] = interface_times[i];
        new_params[3 + i] = params[i];
      }
      interface.update_block_params(block_params.first, new_params);
    }
    std::cout << "[ok  ] update_block_params\n"; flush_now();

    // Set up vectors to run svZeroD simulation
    std::vector<double> solutions(interface.system_size_ * interface.num_output_steps_);
    std::vector<double> times(interface.num_output_steps_);
    int error_code = 0;

    // Run svZeroD simulation
    std::cout << "[step] update_state...\n"; flush_now();
    interface.update_state(init_state_y, init_state_ydot);
    std::cout << "[ok  ] update_state\n"; flush_now();

    std::cout << "[step] run_simulation (this may take a moment)...\n"; flush_now();
    Heartbeat hb; // prints dots if SVZERO_DEBUG=1 while we’re inside
    interface.run_simulation(0.0, times, solutions, error_code);
    std::cout << "\n[ok  ] run_simulation returned (error_code=" << error_code << ")\n"; flush_now();
    if (error_code != 0) {
      throw std::runtime_error("run_simulation returned non-zero error_code");
    }

    // Parse output and calculate mean flow/pressure in aorta and coronary
    int sol_idx = 0;
    double mean_aortic_flow = 0.0;
    double mean_aortic_pressure = 0.0;
    double mean_bc_lca1_outlet_flow = 0.0;
    double mean_bc_lca1_outlet_pressure = 0.0;
    for (int tstep = 0; tstep < interface.num_output_steps_; tstep++) {
      for (int state = 0; state < interface.system_size_; state++) {
        sol_idx = interface.system_size_ * tstep + state;
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
    }
    mean_aortic_flow /= (double)interface.num_output_steps_;
    mean_aortic_pressure /= (double)interface.num_output_steps_;
    mean_bc_lca1_outlet_flow /= (double)interface.num_output_steps_;
    mean_bc_lca1_outlet_pressure /= (double)interface.num_output_steps_;

    std::cout << "Simulation output:\n";
    std::cout << "Mean aortic flow = " << mean_aortic_flow << "\n";
    std::cout << "Mean aortic pressure = " << mean_aortic_pressure << "\n";
    std::cout << "Mean BC_lca1 outlet flow = " << mean_bc_lca1_outlet_flow << "\n";
    std::cout << "Mean BC_lca1 outlet pressure = " << mean_bc_lca1_outlet_pressure << "\n";
    flush_now();

    // Compare with pre-computed values
    double error_limit = 0.05;
    std::vector<std::string> wrong_quantities;
    bool is_wrong = false;
    auto rel_err = [](double v, double ref){ return std::abs(v / ref - 1.0); };
    if (rel_err(mean_aortic_flow, 268.23) > error_limit) {
      is_wrong = true; wrong_quantities.push_back("Mean aortic flow");
    }
    if (rel_err(mean_aortic_pressure, 113.443) > error_limit) {
      is_wrong = true; wrong_quantities.push_back("Mean aortic pressure");
    }
    if (rel_err(mean_bc_lca1_outlet_flow, 0.00755739) > error_limit) {
      is_wrong = true; wrong_quantities.push_back("Mean BC_lca1 outlet flow");
    }
    if (rel_err(mean_bc_lca1_outlet_pressure, 5.88295) > error_limit) {
      is_wrong = true; wrong_quantities.push_back("Mean BC_lca1 outlet pressure");
    }

    if (is_wrong) {
      std::string error_msg = " ";
      for (const auto& q : wrong_quantities) error_msg += q + "; ";
      throw std::runtime_error("Error in the following quantities:" + error_msg);
    }

    std::cout << "[PASS] interface test 01\n"; flush_now();
    return 0;

  } catch (const std::exception& e) {
    std::cerr << "[EXC] " << e.what() << "\n"; flush_now();
    return 1;
  } catch (...) {
    std::cerr << "[EXC] unknown exception\n"; flush_now();
    return 1;
  }
}
