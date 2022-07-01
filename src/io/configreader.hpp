#ifndef SVZERODSOLVER_IO_CONFIGREADER_HPP_
#define SVZERODSOLVER_IO_CONFIGREADER_HPP_

#include <string>
#include <stdexcept>

#include "../model/model.hpp"
#include "../model/junction.hpp"
#include "../helpers/debug.hpp"
#include "../external/simdjson/singleheader/simdjson.h"

namespace IO
{

    template <typename T>
    class ConfigReader
    {
    public:
        ConfigReader();
        ConfigReader(std::string filename);
        ~ConfigReader();
        MODEL::Model<T> get_model();
        int get_num_time_steps();
        T get_time_step_size();
        T cardiac_cycle_period = 1.0;

    private:
        bool model_created = false;
        simdjson::dom::parser parser;
        simdjson::dom::element config;
        static bool has_key(simdjson::dom::element &ele, std::string key);
        static T get_default(simdjson::dom::element &ele, std::string key, T def);
        static simdjson::dom::element get_default(simdjson::dom::element &ele, std::string key, simdjson::dom::element def);
        static MODEL::TimeDependentParameter<T> get_time_dependent_parameter(simdjson::dom::element &json_times, simdjson::dom::element &json_values);
    };

    template <typename T>
    ConfigReader<T>::ConfigReader(std::string filename)
    {
        config = parser.load(filename);
    }

    template <typename T>
    ConfigReader<T>::~ConfigReader()
    {
    }

    template <typename T>
    bool ConfigReader<T>::has_key(simdjson::dom::element &ele, std::string key)
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

    template <typename T>
    T ConfigReader<T>::get_default(simdjson::dom::element &ele, std::string key, T def)
    {
        T value = def;
        if (has_key(ele, key))
        {
            value = ele[key];
        }
        return value;
    }

    template <typename T>
    simdjson::dom::element ConfigReader<T>::get_default(simdjson::dom::element &ele, std::string key, simdjson::dom::element def)
    {
        simdjson::dom::element value = def;
        if (has_key(ele, key))
        {
            value = ele[key];
        }
        return value;
    }

    template <typename T>
    MODEL::TimeDependentParameter<T> ConfigReader<T>::get_time_dependent_parameter(simdjson::dom::element &json_times, simdjson::dom::element &json_values)
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
        return MODEL::TimeDependentParameter<T>(times, values);
    }

    template <typename T>
    MODEL::Model<T> ConfigReader<T>::get_model()
    {
        // Create blog mapping
        MODEL::Model<T> model;

        // Create list to store block connections while generating blocks
        std::vector<std::tuple<std::string, std::string>> connections;

        // Create junctions
        for (simdjson::dom::element junction_config : config["junctions"])
        {
            if ((static_cast<std::string>(junction_config["junction_type"]) == "NORMAL_JUNCTION") || (static_cast<std::string>(junction_config["junction_type"]) == "internal_junction"))
            {
                std::string name = static_cast<std::string>(junction_config["junction_name"]);
                model.blocks.insert(std::make_pair(name, MODEL::Junction<T>(name)));
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
                model.blocks.insert(std::make_pair(name, MODEL::BloodVessel<T>(R = R, C = C, L = L, stenosis_coefficient = stenosis_coefficient, name = name)));
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
                                    model.blocks.insert(std::make_pair(bc_name, MODEL::RCRBlockWithDistalPressure<T>(Rp = Rp, C = C, Rd = Rd, Pd = Pd, bc_name)));
                                    DEBUG_MSG("Created boundary condition " << bc_name);
                                }
                                else if (static_cast<std::string>(bc_config["bc_type"]) == "FLOW")
                                {
                                    simdjson::dom::element Q_json = bc_values["Q"];
                                    simdjson::dom::element t_json = get_default(bc_values, "t", simdjson::dom::element());
                                    auto Q = get_time_dependent_parameter(t_json, Q_json);
                                    cardiac_cycle_period = Q.cycle_period;
                                    model.blocks.insert(std::make_pair(bc_name, MODEL::FlowReference<T>(Q = Q, bc_name)));
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
                                    model.blocks.insert(std::make_pair(bc_name, MODEL::RCRBlockWithDistalPressure<T>(Rp = Rp, C = C, Rd = Rd, Pd = Pd, bc_name)));
                                    DEBUG_MSG("Created boundary condition " << bc_name);
                                }
                                else if (static_cast<std::string>(bc_config["bc_type"]) == "FLOW")
                                {
                                    simdjson::dom::element Q_json = bc_values["Q"];
                                    simdjson::dom::element t_json = get_default(bc_values, "t", simdjson::dom::element());
                                    auto Q = get_time_dependent_parameter(t_json, Q_json);
                                    cardiac_cycle_period = Q.cycle_period;
                                    model.blocks.insert(std::make_pair(bc_name, MODEL::FlowReference<T>(Q = Q, bc_name)));
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
                                    {if ((ele1.name == std::get<0>(connection)) && (ele2.name == std::get<1>(connection))){ model.nodes.push_back(MODEL::Node(ele1.name + "_" + ele2.name)); DEBUG_MSG("Created node " << model.nodes.back().name); ele1.outlet_nodes.push_back(&model.nodes.back()); ele2.inlet_nodes.push_back(&model.nodes.back()); model.nodes.back().setup_dofs(model.dofhandler); } },
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
        model_created = true;
        return model;
    }

    template <typename T>
    int ConfigReader<T>::get_num_time_steps()
    {
        T num_cycles = config["simulation_parameters"]["number_of_cardiac_cycles"];
        T num_pts_per_cycle = config["simulation_parameters"]["number_of_time_pts_per_cardiac_cycle"];
        int num_time_steps = (num_pts_per_cycle - 1) * num_cycles + 1;
        return num_time_steps;
    }

    template <typename T>
    T ConfigReader<T>::get_time_step_size()
    {
        if (model_created == false)
        {
            // Cardiac cycle period is set at model creation so it has to be performed first
            throw std::runtime_error("Please create model before calculating timstepsize.");
        }
        T num_pts_per_cycle = config["simulation_parameters"]["number_of_time_pts_per_cardiac_cycle"];
        T time_step_size = cardiac_cycle_period / (num_pts_per_cycle - 1);
        return time_step_size;
    }
} // namespace IO

#endif // SVZERODSOLVER_IO_CONFIGREADER_HPP_