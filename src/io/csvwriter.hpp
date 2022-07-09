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
 * @param path Path to the output csv file
 * @param times Sequence of time steps corresponding to the solutions
 * @param states Sequence of states corresponding to the time steps
 * @param model The underlying model
 * @param mean Toggle whether only the mean over all time steps should be
 * written
 */
template <typename T>
void write_csv(std::string &path, std::vector<T> &times,
               std::vector<ALGEBRA::State<T>> &states, MODEL::Model<T> &model,
               bool mean = false) {
  // Create string stream to buffer output
  std::stringstream out;

  // Write column labels
  out << "name,time,flow_in,flow_out,pressure_in,pressure_out\n";

  // Create buffer for lines
  char buff[100];

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
        sprintf(buff, "%s,,%.10e,%.10e,%.10e,%.10e\n", name.c_str(),
                inflow_mean, outflow_mean, inpres_mean, outpres_mean);
        out << buff;
      } else {
        for (size_t i = 0; i < times.size(); i++) {
          sprintf(buff, "%s,%.10f,%.10e,%.10e,%.10e,%.10e\n", name.c_str(),
                  times[i], states[i].y[inflow_dof], states[i].y[outflow_dof],
                  states[i].y[inpres_dof], states[i].y[outpres_dof]);
          out << buff;
        }
      }
    }
  }

  std::ofstream ofs(path);
  ofs << out.str();
  ofs.close();
}
}  // namespace IO

#endif  // SVZERODSOLVER_IO_CSVWRITER_HPP_