// #include <json/value.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <list>
#include <variant>

#include "junction.hpp"
#include "bloodvessel.hpp"
#include "rcrblockwithdistalpressure.hpp"
#include "flowreference.hpp"
#include "node.hpp"
#include "integrator.hpp"
#include "model.hpp"
#include "system.hpp"
#include "parameter.hpp"
#include "reader.hpp"
#include "writer.hpp"

#include <stdexcept>
#ifdef DEBUG
#define DEBUG_MSG(str)                 \
    do                                 \
    {                                  \
        std::cout << str << std::endl; \
    } while (false)
#else
#define DEBUG_MSG(str) \
    do                 \
    {                  \
    } while (false)
#endif

typedef double T;

template <typename TT>
using S = SparseSystem<TT>;

int main(int argc, char *argv[])
{

    std::string input_file = argv[1];
    std::string output_file = argv[2];
    bool mean = false;
    if (argc > 3)
    {
        mean = true;
    }
    bool steady_inital = true;
    int max_iter = 30;
    int output_interval = 1;

    DEBUG_MSG("Starting svZeroDSolver");
    DEBUG_MSG("Reading configuration from " << input_file);
    ConfigReader<T> config(input_file);

    // Create the blocks and get simulation parameters
    DEBUG_MSG("Creating model");
    auto model = config.get_model();
    T time_step_size = config.get_time_step_size();
    int num_time_steps = config.get_num_time_steps();

    DEBUG_MSG("Setup simulutation");
    DEBUG_MSG("Number of timesteps: " << num_time_steps);
    DEBUG_MSG("Time step size:      " << time_step_size);
    DEBUG_MSG("Size of system:      " << model.dofhandler.size());

    DEBUG_MSG("Starting simulation");
    State<T> state = State<T>::Zero(model.dofhandler.size());
    S<T> system(model.dofhandler.size());
    system.reserve(model.get_num_triplets());

    // Create steady initial
    if (steady_inital)
    {
        DEBUG_MSG("Calculating steady initial condition");
        T time_step_size_steady = config.cardiac_cycle_period / 10.0;
        auto model_steady = config.get_model();
        model_steady.to_steady();
        model_steady.update_constant(system);
        Integrator<T, S> integrator_steady(system, time_step_size_steady, 0.1);
        for (size_t i = 0; i < 31; i++)
        {
            state = integrator_steady.step(state, time_step_size_steady * T(i), model_steady, max_iter);
        }
    }
    model.update_constant(system);

    Integrator<T, S> integrator(system, time_step_size, 0.1);

    std::vector<State<T>> states;
    std::vector<T> times;
    states.reserve(num_time_steps + 1);
    times.reserve(num_time_steps + 1);

    T time = 0.0;

    states.push_back(state);
    times.push_back(time);

    int interval_counter = 0;
    for (size_t i = 0; i < num_time_steps; i++)
    {
        state = integrator.step(state, time, model, max_iter);
        interval_counter += 1;
        time = time_step_size * T(i + 1);
        if (interval_counter == output_interval)
        {
            times.push_back(time);
            states.push_back(state);
            interval_counter = 0;
        }
    }
    DEBUG_MSG("Simulation completed");

    write_csv<T>(output_file, times, states, model, mean);

    return 0;
}
