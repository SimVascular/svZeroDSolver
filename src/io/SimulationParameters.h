#ifndef SVZERODSOLVER_SIMULATIONPARAMETERS_HPP_
#define SVZERODSOLVER_SIMULATIONPARAMETERS_HPP_

#include <list>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>

//#include "../helpers/debug.hpp"
#include "../helpers/helpers.h"
#include "../model/Model.h"
#include "../algebra/State.h"

namespace io {

/**
 * @brief Simulation parameters
 *
 * @tparam T Scalar type (e.g. `float`, `double`)
 */
struct SimulationParameters 
{
  // Negative value indicates this has not
  // been read from config file yet.
  double sim_time_step_size;  ///< Simulation time step size
  double sim_abs_tol;         ///< Absolute tolerance for simulation

  int sim_num_cycles;      ///< Number of cardiac cycles to simulate
  int sim_pts_per_cycle;   ///< Number of time steps per cardiac cycle
  int sim_num_time_steps;  ///< Total number of time steps
  int sim_nliter;          ///< Maximum number of non-linear iterations in time
                           ///< integration
  int output_interval;     ///< Interval of writing output

  bool sim_steady_initial;  ///< Start from steady solution
  bool output_variable_based;  ///< Output variable based instead of vessel based
  bool output_mean_only;      ///< Output only the mean value
  bool output_derivative;     ///< Output derivatives
  bool output_all_cycles;     ///< Output all cardiac cycles

  bool sim_coupled;  ///< Running 0D simulation coupled with external solver
  double sim_external_step_size;  ///< Step size of external solver if running
                             ///< coupled
};

algebra::State load_initial_condition(const nlohmann::json& config, zd_model::Model& model);

void load_simulation_model(const nlohmann::json& config, zd_model::Model& model);

SimulationParameters load_simulation_params(const nlohmann::json& config);


}  

#endif  
