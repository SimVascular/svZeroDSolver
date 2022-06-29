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
#include "writer.hpp"
#include "simdjson.h"

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
using S = System<TT>;

TimeDependentParameter<T> get_time_dependent_parameter(simdjson::dom::element &json_times, simdjson::dom::element &json_values)
{
    std::vector<T> times;
    std::vector<T> values;
    if (json_values.is_double())
    {
        values.push_back(json_values);
    }
    else
    {
        int len = simdjson::dom::array(json_times).size();
        times.reserve(len);
        values.reserve(len);
        for (T time : json_times)
        {
            times.push_back(time);
        }
        for (T value : json_values)
        {
            values.push_back(value);
        }
    }
    return TimeDependentParameter<T>(times, values);
}

bool has_key(simdjson::dom::element &ele, std::string key)
{
    bool has_key = true;
    try
    {
        simdjson::dom::element tmp = ele.at_key(key);
    }
    catch (simdjson::simdjson_error)
    {
        has_key = false;
    }
    return has_key;
}

T get_default(simdjson::dom::element &ele, std::string key, T def)
{
    T value = def;
    if (has_key(ele, key))
    {
        value = ele[key];
    }
    return value;
}

simdjson::dom::element get_default(simdjson::dom::element &ele, std::string key, simdjson::dom::element def)
{
    simdjson::dom::element value = def;
    if (has_key(ele, key))
    {
        value = ele[key];
    }
    return value;
}

Model<T> create_model(simdjson::dom::element &config)
{
    // Create blog mapping
    Model<T> model;

    // Create list to store block connections while generating blocks
    std::vector<std::tuple<std::string, std::string>> connections;

    // Create junctions
    for (simdjson::dom::element junction_config : config["junctions"])
    {
        if ((static_cast<std::string>(junction_config["junction_type"]) == "NORMAL_JUNCTION") || (static_cast<std::string>(junction_config["junction_type"]) == "internal_junction"))
        {
            std::string name = static_cast<std::string>(junction_config["junction_name"]);
            model.blocks.insert(std::make_pair(name, Junction<T>(name)));
            DEBUG_MSG("Created junction " << name);
            for (simdjson::dom::element inlet_vessel : junction_config["inlet_vessels"])
            {
                auto connection = std::make_tuple("V" + std::to_string(inlet_vessel.get_int64()), name);
                connections.push_back(connection);
                DEBUG_MSG("Found connection " << std::get<0>(connection) << "/" << std::get<1>(connection));
            }
            for (simdjson::dom::element outlet_vessel : junction_config["outlet_vessels"])
            {
                auto connection = std::make_tuple(name, "V" + std::to_string(outlet_vessel.get_int64()));
                connections.push_back(connection);
                DEBUG_MSG("Found connection " << std::get<0>(connection) << "/" << std::get<1>(connection));
            }
        }
        else
        {
            throw std::invalid_argument("Unknown junction type " + static_cast<std::string>(junction_config["junction_type"]));
        }
    }

    // Create vessels
    for (simdjson::dom::element vessel_config : config["vessels"])
    {
        if (static_cast<std::string>(vessel_config["zero_d_element_type"]) == "BloodVessel")
        {
            simdjson::dom::element vessel_values = vessel_config["zero_d_element_values"];
            std::string vessel_id_string = std::to_string(vessel_config["vessel_id"].get_int64());
            std::string name = "V" + vessel_id_string;
            T R = vessel_values["R_poiseuille"];
            T C = get_default(vessel_values, "C", 0.0);
            T L = get_default(vessel_values, "L", 0.0);
            T stenosis_coefficient = get_default(vessel_values, "stenosis_coefficient", 0.0);
            model.blocks.insert(std::make_pair(name, BloodVessel<T>(R = R, C = C, L = L, stenosis_coefficient = stenosis_coefficient, name = name)));
            DEBUG_MSG("Created vessel " << name);

            if (has_key(vessel_config, "boundary_conditions") == true)
            {
                simdjson::dom::element vessel_bcs = vessel_config["boundary_conditions"];
                if (has_key(vessel_bcs, "inlet"))
                {
                    std::string bc_name = "BC" + vessel_id_string + "_inlet";
                    std::string bc_handle = static_cast<std::string>(vessel_bcs["inlet"]);
                    for (simdjson::dom::element bc_config : config["boundary_conditions"])
                    {
                        if (static_cast<std::string>(bc_config["bc_name"]) == bc_handle)
                        {
                            simdjson::dom::element bc_values = bc_config["bc_values"];
                            if (static_cast<std::string>(bc_config["bc_type"]) == "RCR")
                            {
                                T Rp = bc_values["Rp"];
                                T C = bc_values["C"];
                                T Rd = bc_values["Rd"];
                                T Pd = bc_values["Pd"];
                                model.blocks.insert(std::make_pair(bc_name, RCRBlockWithDistalPressure<T>(Rp = Rp, C = C, Rd = Rd, Pd = Pd, bc_name)));
                                DEBUG_MSG("Created boundary condition " << bc_name);
                            }
                            else if (static_cast<std::string>(bc_config["bc_type"]) == "FLOW")
                            {
                                simdjson::dom::element Q_json = bc_values["Q"];
                                simdjson::dom::element t_json = get_default(bc_values, "t", simdjson::dom::element());
                                auto Q = get_time_dependent_parameter(t_json, Q_json);
                                model.time_params.push_back(Q);
                                model.blocks.insert(std::make_pair(bc_name, FlowReference<T>(Q = Q, bc_name)));
                                DEBUG_MSG("Created boundary condition " << bc_name);
                            }
                            else
                            {
                                throw std::invalid_argument("Unknown boundary condition type " + static_cast<std::string>(bc_config["bc_type"]));
                            }
                            break;
                        }
                    }
                    auto connection = std::make_tuple(bc_name, name);
                    connections.push_back(connection);
                    DEBUG_MSG("Found connection " << std::get<0>(connection) << "/" << std::get<1>(connection));
                }
                if (has_key(vessel_bcs, "outlet") == true)
                {
                    std::string bc_name = "BC" + vessel_id_string + "_outlet";
                    std::string bc_handle = static_cast<std::string>(vessel_bcs["outlet"]);
                    for (simdjson::dom::element bc_config : config["boundary_conditions"])
                    {
                        if (static_cast<std::string>(bc_config["bc_name"]) == bc_handle)
                        {
                            simdjson::dom::element bc_values = bc_config["bc_values"];
                            if (static_cast<std::string>(bc_config["bc_type"]) == "RCR")
                            {
                                T Rp = bc_values["Rp"];
                                T C = bc_values["C"];
                                T Rd = bc_values["Rd"];
                                T Pd = bc_values["Pd"];
                                model.blocks.insert(std::make_pair(bc_name, RCRBlockWithDistalPressure<T>(Rp = Rp, C = C, Rd = Rd, Pd = Pd, bc_name)));
                                DEBUG_MSG("Created boundary condition " << bc_name);
                            }
                            else if (static_cast<std::string>(bc_config["bc_type"]) == "FLOW")
                            {
                                simdjson::dom::element Q_json = bc_values["Q"];
                                simdjson::dom::element t_json = get_default(bc_values, "t", simdjson::dom::element());
                                auto Q = get_time_dependent_parameter(t_json, Q_json);
                                model.time_params.push_back(Q);
                                model.blocks.insert(std::make_pair(bc_name, FlowReference<T>(Q = Q, bc_name)));
                                DEBUG_MSG("Created boundary condition " << bc_name);
                            }
                            else
                            {
                                throw std::invalid_argument("Unknown boundary condition type " + static_cast<std::string>(bc_config["bc_type"]));
                            }
                            break;
                        }
                    }
                    auto connection = std::make_tuple(name, bc_name);
                    connections.push_back(connection);
                    DEBUG_MSG("Found connection " << std::get<0>(connection) << "/" << std::get<1>(connection));
                }
            }
        }
        else
        {
            throw std::invalid_argument("Unknown vessel type " + static_cast<std::string>(vessel_config["vessel_type"]));
        }
    }

    // Create Connections
    for (auto &connection : connections)
    {
        for (auto &[key, elem1] : model.blocks)
        {
            std::visit([&](auto &&ele1)
                       {
                        for (auto &[key, elem2] : model.blocks)
                        {
                            std::visit([&](auto &&ele2)
                                    {if ((ele1.name == std::get<0>(connection)) && (ele2.name == std::get<1>(connection))){ model.nodes.push_back(Node(ele1.name + "_" + ele2.name)); DEBUG_MSG("Created node " << model.nodes.back().name); ele1.outlet_nodes.push_back(&model.nodes.back()); ele2.inlet_nodes.push_back(&model.nodes.back()); model.nodes.back().setup_dofs(model.dofhandler); } },
                                    elem2);
                        } },
                       elem1);
        }
    }
    for (auto &[key, elem] : model.blocks)
    {
        std::visit([&](auto &&block)
                   { block.setup_dofs(model.dofhandler); },
                   elem);
    }
    return model;
}

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
    simdjson::dom::parser parser;
    simdjson::dom::element config = parser.load(input_file);

    // Create the blocks
    DEBUG_MSG("Creating model");
    auto model = create_model(config);

    simdjson::dom::element sim_params = config["simulation_parameters"];
    T cardiac_cycle_period = 1.0;
    for (auto tparam : model.time_params)
    {
        cardiac_cycle_period = tparam.cycle_period;
        break;
    }
    T num_cycles = sim_params["number_of_cardiac_cycles"];
    T num_pts_per_cycle = sim_params["number_of_time_pts_per_cardiac_cycle"];
    T time_step_size = cardiac_cycle_period / (num_pts_per_cycle - 1);
    int num_time_steps = (num_pts_per_cycle - 1) * num_cycles + 1;

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
        T time_step_size_steady = cardiac_cycle_period / 10.0;
        auto model_steady = create_model(config);
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
