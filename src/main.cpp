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

TimeDependentParameter<T> get_time_dependent_parameter(Json::Value &json_times, Json::Value &json_values)
{
    std::vector<T> times;
    std::vector<T> values;
    if (json_values.isDouble())
    {
        values.push_back(json_values.asDouble());
    }
    else
    {
        for (auto it = json_times.begin(); it != json_times.end(); it++)
        {
            times.push_back(it->asDouble());
        }
        for (auto it = json_values.begin(); it != json_values.end(); it++)
        {
            values.push_back(it->asDouble());
        }
    }
    return TimeDependentParameter<T>(times, values);
}

Model<T> create_model(Json::Value &config)
{
    // Create blog mapping
    Model<T> model;

    // Create list to store block connections while generating blocks
    std::vector<std::tuple<std::string, std::string>> connections;

    // Create junctions
    int num_junctions = config["junctions"].size();
    DEBUG_MSG("Number of junctions: " << num_junctions);
    for (int i = 0; i < num_junctions; i++)
    {
        if ((config["junctions"][i]["junction_type"].asString() == "NORMAL_JUNCTION") || (config["junctions"][i]["junction_type"].asString() == "internal_junction"))
        {
            std::string name = config["junctions"][i]["junction_name"].asString();
            model.blocks.insert(std::make_pair(name, Junction<T>(name)));
            DEBUG_MSG("Created junction " << name);
            for (int j = 0; j < config["junctions"][i]["inlet_vessels"].size(); j++)
            {
                auto connection = std::make_tuple("V" + config["junctions"][i]["inlet_vessels"][j].asString(), name);
                connections.push_back(connection);
                DEBUG_MSG("Found connection " << std::get<0>(connection) << "/" << std::get<1>(connection));
            }
            for (int j = 0; j < config["junctions"][i]["outlet_vessels"].size(); j++)
            {
                auto connection = std::make_tuple(name, "V" + config["junctions"][i]["outlet_vessels"][j].asString());
                connections.push_back(connection);
                DEBUG_MSG("Found connection " << std::get<0>(connection) << "/" << std::get<1>(connection));
            }
        }
        else
        {
            throw std::invalid_argument("Unknown junction type " + config["junctions"][i]["junction_type"].asString());
        }
    }

    // Create vessels
    int num_vessels = config["vessels"].size();
    DEBUG_MSG("Number of vessels: " << num_vessels);
    for (int i = 0; i < num_vessels; i++)
    {
        if ((config["vessels"][i]["zero_d_element_type"].asString() == "BloodVessel"))
        {
            Json::Value vessel_values = config["vessels"][i]["zero_d_element_values"];
            std::string name = "V" + config["vessels"][i]["vessel_id"].asString();
            T R = vessel_values["R_poiseuille"].asDouble();
            T C = vessel_values.get("C", 0.0).asDouble();
            T L = vessel_values.get("L", 0.0).asDouble();
            T stenosis_coefficient = vessel_values.get("stenosis_coefficient", 0.0).asDouble();
            model.blocks.insert(std::make_pair(name, BloodVessel<T>(R = R, C = C, L = L, stenosis_coefficient = stenosis_coefficient, name = name)));
            DEBUG_MSG("Created vessel " << name);

            if (config["vessels"][i].isMember("boundary_conditions"))
            {
                if (config["vessels"][i]["boundary_conditions"].isMember("inlet"))
                {
                    std::string bc_name = "BC" + config["vessels"][i]["vessel_id"].asString() + "_inlet";
                    std::string bc_handle = config["vessels"][i]["boundary_conditions"]["inlet"].asString();
                    for (int j = 0; j < config["boundary_conditions"].size(); j++)
                    {
                        if (config["boundary_conditions"][j]["bc_name"].asString() == bc_handle)
                        {
                            Json::Value bc_values = config["boundary_conditions"][j]["bc_values"];
                            if (config["boundary_conditions"][j]["bc_type"].asString() == "RCR")
                            {
                                T Rp = bc_values["Rp"].asDouble();
                                T C = bc_values["C"].asDouble();
                                T Rd = bc_values["Rd"].asDouble();
                                T Pd = bc_values["Pd"].asDouble();
                                model.blocks.insert(std::make_pair(bc_name, RCRBlockWithDistalPressure<T>(Rp = Rp, C = C, Rd = Rd, Pd = Pd, bc_name)));
                                DEBUG_MSG("Created boundary condition " << bc_name);
                            }
                            else if (config["boundary_conditions"][j]["bc_type"].asString() == "FLOW")
                            {
                                Json::Value Q_json = bc_values["Q"];
                                Json::Value t_json = bc_values.get("t", 0.0);
                                auto Q = get_time_dependent_parameter(t_json, Q_json);
                                if (Q.isconstant == false)
                                {
                                    config["simulation_parameters"]["cardiac_cycle_period"] = Json::Value(Q.cycle_period);
                                }
                                model.blocks.insert(std::make_pair(bc_name, FlowReference<T>(Q = Q, bc_name)));
                                DEBUG_MSG("Created boundary condition " << bc_name);
                            }
                            else
                            {
                                throw std::invalid_argument("Unknown boundary condition type " + config["boundary_conditions"][j]["bc_type"].asString());
                            }
                            break;
                        }
                    }
                    auto connection = std::make_tuple(bc_name, name);
                    connections.push_back(connection);
                    DEBUG_MSG("Found connection " << std::get<0>(connection) << "/" << std::get<1>(connection));
                }
                if (config["vessels"][i]["boundary_conditions"].isMember("outlet"))
                {
                    std::string bc_name = "BC" + config["vessels"][i]["vessel_id"].asString() + "_outlet";
                    std::string bc_handle = config["vessels"][i]["boundary_conditions"]["outlet"].asString();
                    for (int j = 0; j < config["boundary_conditions"].size(); j++)
                    {
                        if (config["boundary_conditions"][j]["bc_name"].asString() == bc_handle)
                        {
                            Json::Value bc_values = config["boundary_conditions"][j]["bc_values"];
                            if (config["boundary_conditions"][j]["bc_type"].asString() == "RCR")
                            {
                                T Rp = bc_values["Rp"].asDouble();
                                T C = bc_values["C"].asDouble();
                                T Rd = bc_values["Rd"].asDouble();
                                T Pd = bc_values["Pd"].asDouble();
                                model.blocks.insert(std::make_pair(bc_name, RCRBlockWithDistalPressure<T>(Rp = Rp, C = C, Rd = Rd, Pd = Pd, bc_name)));
                                DEBUG_MSG("Created boundary condition " << bc_name);
                            }
                            else if (config["boundary_conditions"][j]["bc_type"].asString() == "FLOW")
                            {
                                Json::Value Q_json = bc_values["Q"];
                                Json::Value t_json = bc_values.get("t", 0.0);
                                auto Q = get_time_dependent_parameter(t_json, Q_json);
                                if (Q.isconstant == false)
                                {
                                    config["simulation_parameters"]["cardiac_cycle_period"] = Json::Value(Q.cycle_period);
                                }
                                model.blocks.insert(std::make_pair(bc_name, FlowReference<T>(Q = Q, bc_name)));
                                DEBUG_MSG("Created boundary condition " << bc_name);
                            }
                            else
                            {
                                throw std::invalid_argument("Unknown boundary condition type " + config["boundary_conditions"][j]["bc_type"].asString());
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
            throw std::invalid_argument("Unknown vessel type " + config["vessels"][i]["vessel_type"].asString());
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
    bool steady_inital = true;
    int max_iter = 30;
    int output_interval = 1;

    DEBUG_MSG("Starting svZeroDSolver");
    DEBUG_MSG("Reading configuration from " << input_file);
    std::ifstream file_input(input_file);
    Json::Reader reader;
    Json::Value config;
    reader.parse(file_input, config);

    // Create the blocks
    DEBUG_MSG("Creating model");
    auto model = create_model(config);

    Json::Value sim_params = config["simulation_parameters"];
    T cardiac_cycle_period = sim_params.get("cardiac_cycle_period", 1.0).asDouble();
    T num_cycles = sim_params["number_of_cardiac_cycles"].asDouble();
    T num_pts_per_cycle = sim_params["number_of_time_pts_per_cardiac_cycle"].asDouble();
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

    write_csv<T>(output_file, times, states, model);

    return 0;
}
