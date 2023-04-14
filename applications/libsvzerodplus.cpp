// Copyright (c) Stanford University, The Regents of the University of
//               California, and others.
//
// All Rights Reserved.
//
// See Copyright-SimVascular.txt for additional details.
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject
// to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
// OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/**
 * @file python.cpp
 * @brief Python binding for svZeroDSolver
 */
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <io/jsonhandler.hpp>

#include "solve/solver.hpp"
#include "pybind11_json/pybind11_json.hpp"

namespace py = pybind11;

PYBIND11_MODULE(libsvzerodplus, m) {
  using Solver = SOLVE::Solver<double>;
  m.doc() = "svZeroDSolver";
  py::class_<Solver>(m, "Solver")
    .def(py::init([] (py::dict& config) {
      nlohmann::json config_json = config;
      auto config_handler = IO::JsonHandler(config_json);
      return new Solver(config_handler);
      }))
    .def("run",&Solver::run)
    .def("get_full_result",&Solver::get_full_result)
    .def("get_single_result",&Solver::get_single_result)
    .def("get_single_result_avg",&Solver::get_single_result_avg)
    .def("update_block_params",&Solver::update_block_params);
}
