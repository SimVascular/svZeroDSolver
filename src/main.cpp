// #include <json/value.h>
#include <fstream>
#include <iostream>
#include "json.h"
#include <vector>
#include <map>
#include <list>
#include "Eigen/Dense"

#include "block.hpp"
#include "junction.hpp"
#include "bloodvessel.hpp"

#include <stdexcept>

std::map<std::string, Block *> create_blocks(Json::Value &config)
{
    // Create blog mapping
    std::map<std::string, Block *> blocks;

    // Create list to store block connections while generating blocks
    std::list<std::tuple<std::string, std::string>> connections;

    // Create junctions
    int num_junctions = config["junctions"].size();
    std::cout << "Number of junctions: " << num_junctions << std::endl;
    for (int i = 0; i < num_junctions; i++)
    {
        if ((config["junctions"][i]["junction_type"].asString() == "NORMAL_JUNCTION") || (config["junctions"][i]["junction_type"].asString() == "internal_junction"))
        {
            Junction::Parameters params;
            Junction junction = Junction(params, config["junctions"][i]["junction_name"].asString());
            blocks[junction.name] = &junction;
            std::cout << "Created junction " << junction.name << std::endl;
            for (int j = 0; j < config["junctions"][i]["inlet_vessels"].size(); j++)
            {
                auto connection = std::make_tuple("V" + config["junctions"][i]["inlet_vessels"][j].asString(), junction.name);
                connections.push_back(connection);
                std::cout << "Found connection " << std::get<0>(connection) << "/" << std::get<1>(connection) << std::endl;
            }
            for (int j = 0; j < config["junctions"][i]["outlet_vessels"].size(); j++)
            {
                auto connection = std::make_tuple(junction.name, "V" + config["junctions"][i]["outlet_vessels"][j].asString());
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
            BloodVessel::Parameters params;
            BloodVessel vessel = BloodVessel(params, "V" + config["vessels"][i]["vessel_id"].asString());
            blocks[vessel.name] = &vessel;
            std::cout << "Created vessel " << vessel.name << std::endl;

            if (config["vessels"][i].isMember("boundary_conditions"))
            {
                if (config["vessels"][i]["boundary_conditions"].isMember("inlet"))
                {
                    auto connection = std::make_tuple("BC" + config["vessels"][i]["vessel_id"].asString() + "_inlet", vessel.name);
                    connections.push_back(connection);
                    std::cout << "Found connection " << std::get<0>(connection) << "/" << std::get<1>(connection) << std::endl;
                }
                if (config["vessels"][i]["boundary_conditions"].isMember("outlet"))
                {
                    auto connection = std::make_tuple("BC" + config["vessels"][i]["vessel_id"].asString() + "_outlet", vessel.name);
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

    // Create Boundary Conditions

    return blocks;
}

int main(int argc, char *argv[])
{
    std::cout << "Starting svZeroDSolver" << std::endl;
    std::cout << "Reading configuration from " << argv[1] << std::endl;
    std::ifstream file_input(argv[1]);
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
    auto blocks = create_blocks(config);
}
