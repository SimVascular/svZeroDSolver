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
 * @file jsonwriter.hpp
 * @brief IO::write_json source file
 */
#ifndef SVZERODSOLVER_IO_JSONWRITER_HPP_
#define SVZERODSOLVER_IO_JSONWRITER_HPP_

#include <json.h>

#include <vector>

#include "../algebra/state.hpp"
#include "../helpers/startswith.hpp"
#include "../model/model.hpp"

namespace IO {

/**
 * @brief Write the solution to a json file
 *
 * @tparam T Scalar type (e.g. `float`, `double`)
 * @param times Sequence of time steps corresponding to the solutions
 * @param states Sequence of states corresponding to the time steps
 * @param model The underlying model
 * @return Json encoded output string
 */
template <typename T>
std::string write_json(std::vector<T> &times,
                       std::vector<ALGEBRA::State<T>> &states,
                       MODEL::Model<T> &model) {
  Json::Value output;
  Json::Value json_times(Json::arrayValue);
  for (auto time : times) {
    json_times.append(Json::Value(time));
  }

  Json::Value json_names(Json::arrayValue);
  Json::Value json_flow_in(Json::arrayValue);
  Json::Value json_flow_out(Json::arrayValue);
  Json::Value json_pres_in(Json::arrayValue);
  Json::Value json_pres_out(Json::arrayValue);

  for (auto &[key, elem] : model.blocks) {
    std::string name = "NoName";
    unsigned int inflow_dof;
    unsigned int outflow_dof;
    unsigned int inpres_dof;
    unsigned int outpres_dof;
    std::visit(
        [&](auto &&block) {
          if (HELPERS::startswith(block.name, "V")) {
            name = block.name;
            inflow_dof = block.inlet_nodes[0]->flow_dof;
            outflow_dof = block.outlet_nodes[0]->flow_dof;
            inpres_dof = block.inlet_nodes[0]->pres_dof;
            outpres_dof = block.outlet_nodes[0]->pres_dof;
          }
        },
        elem);

    if (name != "NoName") {
      json_names.append(name);
      Json::Value json_flow_in_i(Json::arrayValue);
      Json::Value json_flow_out_i(Json::arrayValue);
      Json::Value json_pres_in_i(Json::arrayValue);
      Json::Value json_pres_out_i(Json::arrayValue);
      for (auto state : states) {
        json_flow_in_i.append(state.y[inflow_dof]);
        json_flow_out_i.append(state.y[outflow_dof]);
        json_pres_in_i.append(state.y[inpres_dof]);
        json_pres_out_i.append(state.y[outpres_dof]);
      }
      json_flow_in.append(json_flow_in_i);
      json_flow_out.append(json_flow_out_i);
      json_pres_in.append(json_pres_in_i);
      json_pres_out.append(json_pres_out_i);
    }
  }

  output["time"] = json_times;
  output["names"] = json_names;
  output["flow_in"] = json_flow_in;
  output["flow_out"] = json_flow_out;
  output["pressure_in"] = json_pres_in;
  output["pressure_out"] = json_pres_out;

  Json::FastWriter writer;
  return writer.write(output);
}
}  // namespace IO

#endif  // SVZERODSOLVER_IO_JSONWRITER_HPP_