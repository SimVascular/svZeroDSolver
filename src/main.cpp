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
            // Junction junction = Junction()
            Junction::Parameters params;
            Junction junction = Junction(params, config["junctions"][i]["junction_name"].asString());
            std::cout << "Created junction " << junction.name << std::endl;
        }
    }

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
