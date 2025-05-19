// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause
/**
 * @file pysvzerod.cpp
 * @brief Python interface for svZeroDSolver
 */
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Solver.h"
#include "calibrate.h"
#include "pybind11_json/pybind11_json.hpp"

namespace py = pybind11;

PYBIND11_MODULE(pysvzerod, m) {
  using Solver = Solver;
  py::class_<Solver>(m, "Solver")
      .def(py::init([](py::dict& config) {
        const nlohmann::json& config_json = config;
        return Solver(config_json);
      }))
      .def(py::init([](std::string config_file) {
        std::ifstream ifs(config_file);
        const auto& config_json = nlohmann::json::parse(ifs);
        return Solver(config_json);
      }))
      .def("run", &Solver::run)
      .def("get_times", &Solver::get_times)
      .def("get_single_result", &Solver::get_single_result)
      .def("get_single_result_avg", &Solver::get_single_result_avg)
      .def("update_block_params", &Solver::update_block_params)
      .def("read_block_params", &Solver::read_block_params)
      .def("get_full_result", [](Solver& solver) {
        py::module_ pd = py::module_::import("pandas");
        py::module_ io = py::module_::import("io");
        auto result = solver.get_full_result();
        return pd.attr("read_csv")(io.attr("StringIO")(result));
      });

  m.def("simulate", [](py::dict& config) {
    py::module_ pd = py::module_::import("pandas");
    py::module_ io = py::module_::import("io");
    const nlohmann::json& config_json = config;
    auto solver = Solver(config_json);
    solver.run();
    return pd.attr("read_csv")(io.attr("StringIO")(solver.get_full_result()));
  });
  m.def("simulate", [](std::string config_file) {
    py::module_ pd = py::module_::import("pandas");
    py::module_ io = py::module_::import("io");
    std::ifstream ifs(config_file);
    const auto& config_json = nlohmann::json::parse(ifs);
    auto solver = Solver(config_json);
    solver.run();
    return pd.attr("read_csv")(io.attr("StringIO")(solver.get_full_result()));
  });
  m.def("calibrate", [](py::dict& config) {
    const nlohmann::json& config_json = config;
    return calibrate(config);
  });
  m.def("run_simulation_cli", []() {
    py::module_ sys = py::module_::import("sys");
    auto argv = sys.attr("argv").cast<std::vector<std::string>>();
    if (argv.size() != 3) {
      std::cout
          << "Usage: svzerodsolver path/to/config.json path/to/output.csv"
          << std::endl;
      exit(1);
    }
    std::ifstream ifs(argv[1]);
    const auto& config = nlohmann::json::parse(ifs);
    auto solver = Solver(config);
    solver.run();
    solver.write_result_to_csv(argv[2]);
  });
  m.def("run_calibration_cli", []() {
    py::module_ sys = py::module_::import("sys");
    auto argv = sys.attr("argv").cast<std::vector<std::string>>();
    if (argv.size() != 3) {
      std::cout
          << "Usage: svzerodcalibrator path/to/config.json path/to/output.json"
          << std::endl;
      exit(1);
    }
    std::ifstream ifs(argv[1]);
    const auto& config = nlohmann::json::parse(ifs);
    auto output_config = calibrate(config);
    std::ofstream o(argv[2]);
    o << std::setw(4) << output_config << std::endl;
  });
}
