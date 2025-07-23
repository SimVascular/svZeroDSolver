// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "csv_writer.h"

#include <iomanip>

/**
 * @brief Write results vessel based.
 *
 * @param times Sequence of time steps corresponding to the solutions
 * @param states Sequence of states corresponding to the time steps
 * @param model The underlying model
 * @param mean Toggle whether only the mean over all time steps should be
 * written
 * @param derivative Toggle whether to output time-derivatives
 * @return CSV encoded output string
 */
std::string to_vessel_csv(const std::vector<double> &times,
                          const std::vector<State> &states, const Model &model,
                          bool mean, bool derivative) {
  // Create string stream to buffer output
  std::stringstream out;

  // Set floating point format for the entire stream
  out << std::scientific << std::setprecision(16);

  // Write column labels
  if (derivative) {
    out << "name,time,flow_in,flow_out,pressure_in,pressure_out,d_flow_in,d_"
           "flow_out,d_pressure_in,d_pressure_out\n";
  } else {
    out << "name,time,flow_in,flow_out,pressure_in,pressure_out\n";
  }

  // Determine number of time steps
  int num_steps = times.size();

  int inflow_dof;
  int outflow_dof;
  int inpres_dof;
  int outpres_dof;
  for (size_t i = 0; i < model.get_num_blocks(); i++) {
    auto block = model.get_block(i);
    // Extract global solution indices of the block

    if (dynamic_cast<const BloodVessel *>(block) == nullptr &&
        dynamic_cast<const ChamberSphere *>(block) == nullptr &&
        dynamic_cast<const BloodVesselCRL *>(block) == nullptr &&
        dynamic_cast<const RegazzoniChamber *>(block) == nullptr &&
        dynamic_cast<const RegazzoniValve *>(block) == nullptr) {
      continue;
    }

    std::string name = block->get_name();
    inflow_dof = block->inlet_nodes[0]->flow_dof;
    outflow_dof = block->outlet_nodes[0]->flow_dof;
    inpres_dof = block->inlet_nodes[0]->pres_dof;
    outpres_dof = block->outlet_nodes[0]->pres_dof;

    // Write the solution of the block to the output file
    if (derivative) {
      if (mean) {
        double inflow_mean = 0.0;
        double outflow_mean = 0.0;
        double inpres_mean = 0.0;
        double outpres_mean = 0.0;
        double d_inflow_mean = 0.0;
        double d_outflow_mean = 0.0;
        double d_inpres_mean = 0.0;
        double d_outpres_mean = 0.0;

        for (size_t i = 0; i < num_steps; i++) {
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

        out << name << ",," << inflow_mean << "," << outflow_mean << ","
            << inpres_mean << "," << outpres_mean << "," << d_inflow_mean << ","
            << d_outflow_mean << "," << d_inpres_mean << "," << d_outpres_mean
            << "\n";
      } else {
        for (size_t i = 0; i < num_steps; i++) {
          out << name << "," << times[i] << "," << states[i].y[inflow_dof]
              << "," << states[i].y[outflow_dof] << ","
              << states[i].y[inpres_dof] << "," << states[i].y[outpres_dof]
              << "," << states[i].ydot[inflow_dof] << ","
              << states[i].ydot[outflow_dof] << ","
              << states[i].ydot[inpres_dof] << ","
              << states[i].ydot[outpres_dof] << "\n";
        }
      }
    } else {
      if (mean) {
        double inflow_mean = 0.0;
        double outflow_mean = 0.0;
        double inpres_mean = 0.0;
        double outpres_mean = 0.0;

        for (size_t i = 0; i < num_steps; i++) {
          inflow_mean += states[i].y[inflow_dof];
          outflow_mean += states[i].y[outflow_dof];
          inpres_mean += states[i].y[inpres_dof];
          outpres_mean += states[i].y[outpres_dof];
        }
        inflow_mean /= num_steps;
        outflow_mean /= num_steps;
        inpres_mean /= num_steps;
        outpres_mean /= num_steps;

        out << name << ",," << inflow_mean << "," << outflow_mean << ","
            << inpres_mean << "," << outpres_mean << "\n";
      } else {
        for (size_t i = 0; i < num_steps; i++) {
          out << name << "," << times[i] << "," << states[i].y[inflow_dof]
              << "," << states[i].y[outflow_dof] << ","
              << states[i].y[inpres_dof] << "," << states[i].y[outpres_dof]
              << "\n";
        }
      }
    }
  }

  return out.str();
}

/**
 * @brief Write results variable based.
 *
 * @param times Sequence of time steps corresponding to the solutions
 * @param states Sequence of states corresponding to the time steps
 * @param model The underlying model
 * @param mean Toggle whether only the mean over all time steps should be
 * written
 * @param derivative Toggle whether to output time-derivatives
 * @return CSV encoded output string
 */
std::string to_variable_csv(const std::vector<double> &times,
                            const std::vector<State> &states,
                            const Model &model, bool mean, bool derivative) {
  // Create string stream to buffer output
  std::stringstream out;

  // Set floating point format for the entire stream
  out << std::scientific << std::setprecision(16);

  // Determine number of time steps
  int num_steps = times.size();

  // Write column labels
  if (derivative) {
    out << "name,time,y,ydot\n";
    if (mean) {
      for (size_t i = 0; i < model.dofhandler.size(); i++) {
        std::string name = model.dofhandler.variables[i];
        double mean_y = 0.0;
        double mean_ydot = 0.0;

        for (size_t j = 0; j < num_steps; j++) {
          mean_y += states[j].y[i];
          mean_ydot += states[j].ydot[i];
        }
        mean_y /= num_steps;
        mean_ydot /= num_steps;

        out << name << ",," << mean_y << "," << mean_ydot << "\n";
      }
    } else {
      for (size_t i = 0; i < model.dofhandler.size(); i++) {
        std::string name = model.dofhandler.variables[i];
        for (size_t j = 0; j < num_steps; j++) {
          out << name << "," << times[j] << "," << states[j].y[i] << ","
              << states[j].ydot[i] << "\n";
        }
      }
    }
  } else {
    out << "name,time,y\n";
    if (mean) {
      for (size_t i = 0; i < model.dofhandler.size(); i++) {
        std::string name = model.dofhandler.variables[i];
        double mean_y = 0.0;
        for (size_t j = 0; j < num_steps; j++) {
          mean_y += states[j].y[i];
        }
        mean_y /= num_steps;

        out << name << ",," << mean_y << "\n";
      }
    } else {
      for (size_t i = 0; i < model.dofhandler.size(); i++) {
        std::string name = model.dofhandler.variables[i];
        for (size_t j = 0; j < num_steps; j++) {
          out << name << "," << times[j] << "," << states[j].y[i] << "\n";
        }
      }
    }
  }

  return out.str();
}