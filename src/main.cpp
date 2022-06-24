// #include <json/value.h>
#include <fstream>
#include <iostream>
#include "json.h"
#include <vector>
#include <map>
#include <list>
#include <variant>
#include "Eigen/Dense"

#include "block.hpp"
#include "junction.hpp"
#include "bloodvessel.hpp"
#include "rcrblockwithdistalpressure.hpp"
#include "flowreference.hpp"
#include "node.hpp"
#include "integrator.hpp"
#include "model.hpp"
#include "system.hpp"

#include <stdexcept>

Model create_model(Json::Value &config)
{
    // Create blog mapping
    Model model;

    // Create list to store block connections while generating blocks
    std::vector<std::tuple<std::string, std::string>> connections;

    // Create junctions
    int num_junctions = config["junctions"].size();
    std::cout << "Number of junctions: " << num_junctions << std::endl;
    for (int i = 0; i < num_junctions; i++)
    {
        if ((config["junctions"][i]["junction_type"].asString() == "NORMAL_JUNCTION") || (config["junctions"][i]["junction_type"].asString() == "internal_junction"))
        {
            std::string name = config["junctions"][i]["junction_name"].asString();
            model.blocks.insert(std::make_pair(name, Junction(name)));
            std::cout << "Created junction " << name << std::endl;
            for (int j = 0; j < config["junctions"][i]["inlet_vessels"].size(); j++)
            {
                auto connection = std::make_tuple("V" + config["junctions"][i]["inlet_vessels"][j].asString(), name);
                connections.push_back(connection);
                std::cout << "Found connection " << std::get<0>(connection) << "/" << std::get<1>(connection) << std::endl;
            }
            for (int j = 0; j < config["junctions"][i]["outlet_vessels"].size(); j++)
            {
                auto connection = std::make_tuple(name, "V" + config["junctions"][i]["outlet_vessels"][j].asString());
                connections.push_back(connection);
                std::cout << "Found connection " << std::get<0>(connection) << "/" << std::get<1>(connection) << std::endl;
            }
        }
        else
        {
            throw std::invalid_argument("Unknown junction type " + config["junctions"][i]["junction_type"].asString());
        }
    }

    // Create vessels
    int num_vessels = config["vessels"].size();
    std::cout << "Number of vessels: " << num_vessels << std::endl;
    for (int i = 0; i < num_vessels; i++)
    {
        if ((config["vessels"][i]["zero_d_element_type"].asString() == "BloodVessel"))
        {
            Json::Value vessel_values = config["vessels"][i]["zero_d_element_values"];
            std::string name = "V" + config["vessels"][i]["vessel_id"].asString();
            double R = vessel_values["R_poiseuille"].asDouble();
            double C = vessel_values.get("C", 0.0).asDouble();
            double L = vessel_values.get("L", 0.0).asDouble();
            double stenosis_coefficient = vessel_values.get("stenosis_coefficient", 0.0).asDouble();
            model.blocks.insert(std::make_pair(name, BloodVessel(R = R, C = C, L = L, stenosis_coefficient = stenosis_coefficient, name = name)));
            std::cout << "Created vessel " << name << std::endl;

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
                                double Rp = bc_values["Rp"].asDouble();
                                double C = bc_values["C"].asDouble();
                                double Rd = bc_values["Rd"].asDouble();
                                double Pd = bc_values["Pd"].asDouble();
                                model.blocks.insert(std::make_pair(name, RCRBlockWithDistalPressure(Rp = Rp, C = C, Rd = Rd, Pd = Pd, bc_name)));
                                std::cout << "Created boundary condition " << bc_name << std::endl;
                            }
                            else if (config["boundary_conditions"][j]["bc_type"].asString() == "FLOW")
                            {
                                double Q = bc_values["Q"].asDouble();
                                model.blocks.insert(std::make_pair(bc_name, FlowReference(Q = Q, bc_name)));
                                std::cout << "Created boundary condition " << bc_name << std::endl;
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
                    std::cout << "Found connection " << std::get<0>(connection) << "/" << std::get<1>(connection) << std::endl;
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
                                double Rp = bc_values["Rp"].asDouble();
                                double C = bc_values["C"].asDouble();
                                double Rd = bc_values["Rd"].asDouble();
                                double Pd = bc_values["Pd"].asDouble();
                                model.blocks.insert(std::make_pair(name, RCRBlockWithDistalPressure(Rp = Rp, C = C, Rd = Rd, Pd = Pd, bc_name)));
                                std::cout << "Created boundary condition " << bc_name << std::endl;
                            }
                            else if (config["boundary_conditions"][j]["bc_type"].asString() == "FLOW")
                            {
                                double Q = bc_values["Q"].asDouble();
                                model.blocks.insert(std::make_pair(bc_name, FlowReference(Q = Q, bc_name)));
                                std::cout << "Created boundary condition " << bc_name << std::endl;
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
                    std::cout << "Found connection " << std::get<0>(connection) << "/" << std::get<1>(connection) << std::endl;
                }
            }
        }
        else
        {
            throw std::invalid_argument("Unknown vessel type " + config["vessels"][i]["vessel_type"].asString());
        }
    }

    // Create Connections
    for (size_t i = 0; i < connections.size(); i++)
    {
        auto connection = connections[i];
        for (auto &&elem1 : model.blocks)
        {
            std::visit([&](auto &&ele1)
                       {
                        for (auto &&elem2 : model.blocks)
                        {
                            std::visit([&](auto &&ele2)
                                    {if ((ele1.name == std::get<0>(connection)) && (ele2.name == std::get<1>(connection))){Node node = Node(ele1, ele2, ele1.name + "_" + ele2.name); std::cout << "Created node " << node.name << std::endl; node.setup_dofs(model.dofhandler);}},
                                    elem2.second);
                        } },
                       elem1.second);
        }
    }
    for (auto &&elem : model.blocks)
    {
        std::visit([&](auto &&block)
                   { block.setup_dofs(model.dofhandler); },
                   elem.second);
    }
    return model;
}

int main(int argc, char *argv[])
{
    std::cout << "Starting svZeroDSolver" << std::endl;
    // std::cout << "Reading configuration from " << argv[1] << std::endl;
    std::ifstream file_input("../solver_0d.in");
    Json::Reader reader;
    Json::Value config;
    reader.parse(file_input, config);

    Json::Value sim_params = config["simulation_parameters"];
    double cardiac_cycle_period = sim_params.get("cardiac_cycle_period", 1.0).asDouble();
    double num_cycles = sim_params["number_of_cardiac_cycles"].asDouble();
    double num_pts_per_cycle = sim_params["number_of_time_pts_per_cardiac_cycle"].asDouble();
    double time_step_size = cardiac_cycle_period / (num_pts_per_cycle - 1);
    int num_time_steps = (num_pts_per_cycle - 1) * num_cycles + 1;

    std::cout << "Reading configuration completed" << std::endl;
    std::cout << "Number of timesteps: " << num_time_steps << std::endl;
    std::cout << "Time step size:      " << time_step_size << std::endl;

    // Create the blocks
    auto model = create_model(config);

    System system;
    system.setup_matrices(model.dofhandler.size());

    Integrator integrator = Integrator(system, time_step_size, 0.1);

    State state = State::Zero(model.dofhandler.size());

    double time = 0.0;
    int max_iter = 30;
    integrator.step(state, time, model, max_iter);

    std::cout << system.F << std::endl;

    return 0;
}
