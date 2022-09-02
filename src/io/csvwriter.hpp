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
 * @file csvwriter.hpp
 * @brief IO::write_csv source file
 */
#ifndef SVZERODSOLVER_IO_CSVWRITER_HPP_
#define SVZERODSOLVER_IO_CSVWRITER_HPP_

#include <fstream>
#include <string>
#include <vector>

#include "../algebra/state.hpp"
#include "../helpers/startswith.hpp"
#include "../model/model.hpp"

namespace IO {

/**
 * @brief Write the solution to a csv file
 *
 * @tparam T Scalar type (e.g. `float`, `double`)
 * @param times Sequence of time steps corresponding to the solutions
 * @param states Sequence of states corresponding to the time steps
 * @param model The underlying model
 * @param mean Toggle whether only the mean over all time steps should be
 * written
 * @return CSV encoded output string
 */
template <typename T>
std::string write_csv(std::vector<T> &times,
                      std::vector<ALGEBRA::State<T>> &states,
                      MODEL::Model<T> &model, bool mean = false,
                      bool derivative = false) {
  // Create string stream to buffer output
  std::stringstream out;

  // Create short and long buffer for lines
  char sbuff[130];
  char lbuff[226];

  // Write column labels
  if (derivative) {
    out << "name,time,flow_in,flow_out,pressure_in,pressure_out,d_flow_in,d_"
           "flow_out,d_pressure_in,d_pressure_out\n";
  } else {
    out << "name,time,flow_in,flow_out,pressure_in,pressure_out\n";
  }

  // Determine number of time steps
  T num_steps = times.size();

  unsigned int inflow_dof;
  unsigned int outflow_dof;
  unsigned int inpres_dof;
  unsigned int outpres_dof;
  for (auto &[key, elem] : model.blocks) {
    // Extract global solution indices of the block
    std::string name = "NoName";
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

    // Write the solution of the block to the output file
    if (derivative) {
      if (name != "NoName") {
        if (mean) {
          T inflow_mean = 0.0;
          T outflow_mean = 0.0;
          T inpres_mean = 0.0;
          T outpres_mean = 0.0;
          T d_inflow_mean = 0.0;
          T d_outflow_mean = 0.0;
          T d_inpres_mean = 0.0;
          T d_outpres_mean = 0.0;
          for (size_t i = 0; i < times.size(); i++) {
            inflow_mean += states[i].y[inflow_dof];
            outflow_mean += states[i].y[outflow_dof];
            inpres_mean += states[i].y[inpres_dof];
            outpres_mean += states[i].y[outpres_dof];
            d_inflow_mean += states[i].ydot[inflow_dof];
            d_outflow_mean += states[i].ydot[outflow_dof];
            d_inpres_mean += states[i].ydot[inpres_dof];
            d_outpres_mean += states[i].ydot[outpres_dof];
          }
          inflow_mean /= num_steps;
          outflow_mean /= num_steps;
          inpres_mean /= num_steps;
          outpres_mean /= num_steps;
          d_inflow_mean /= num_steps;
          d_outflow_mean /= num_steps;
          d_inpres_mean /= num_steps;
          d_outpres_mean /= num_steps;
          sprintf(lbuff,
                  "%s,,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n",
                  name.c_str(), inflow_mean, outflow_mean, inpres_mean,
                  outpres_mean, d_inflow_mean, d_outflow_mean, d_inpres_mean,
                  d_outpres_mean);
          out << lbuff;
        } else {
          for (size_t i = 0; i < times.size(); i++) {
            sprintf(
                lbuff,
                "%s,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n",
                name.c_str(), times[i], states[i].y[inflow_dof],
                states[i].y[outflow_dof], states[i].y[inpres_dof],
                states[i].y[outpres_dof], states[i].ydot[inflow_dof],
                states[i].ydot[outflow_dof], states[i].ydot[inpres_dof],
                states[i].ydot[outpres_dof]);
            out << lbuff;
          }
        }
      }
    } else {
      if (name != "NoName") {
        if (mean) {
          T inflow_mean = 0.0;
          T outflow_mean = 0.0;
          T inpres_mean = 0.0;
          T outpres_mean = 0.0;
          for (size_t i = 0; i < times.size(); i++) {
            inflow_mean += states[i].y[inflow_dof];
            outflow_mean += states[i].y[outflow_dof];
            inpres_mean += states[i].y[inpres_dof];
            outpres_mean += states[i].y[outpres_dof];
          }
          inflow_mean /= num_steps;
          outflow_mean /= num_steps;
          inpres_mean /= num_steps;
          outpres_mean /= num_steps;
          sprintf(sbuff, "%s,,%.16e,%.16e,%.16e,%.16e\n", name.c_str(),
                  inflow_mean, outflow_mean, inpres_mean, outpres_mean);
          out << sbuff;
        } else {
          for (size_t i = 0; i < times.size(); i++) {
            sprintf(sbuff, "%s,%.16e,%.16e,%.16e,%.16e,%.16e\n", name.c_str(),
                    times[i], states[i].y[inflow_dof], states[i].y[outflow_dof],
                    states[i].y[inpres_dof], states[i].y[outpres_dof]);
            out << sbuff;
          }
        }
      }
    }
  }

  return out.str();
}
}  // namespace IO

#endif  // SVZERODSOLVER_IO_CSVWRITER_HPP_