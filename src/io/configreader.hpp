/**
 * @file configreader.hpp
 * @brief IO::ConfigReader source file
 */
#ifndef SVZERODSOLVER_IO_CONFIGREADER_HPP_
#define SVZERODSOLVER_IO_CONFIGREADER_HPP_

#include <string>
#include <list>
#include <stdexcept>

#include "../model/model.hpp"
#include "../helpers/debug.hpp"
#include "simdjson.h"

namespace IO
{

    /**
     * @brief svZeroDSolver configuration reader.
     *
     * @tparam T Scalar type (e.g. `float`, `double`)
     */
    template <typename T>
    class ConfigReader
    {
    public:
        /**
         * @brief Construct a new Config Reader object
         */
        ConfigReader();

        /**
         * @brief Construct a new Config Reader object
         *
         * @param filename Name of the configuration file
         */
        ConfigReader(std::string filename);

        /**
         * @brief Destroy the Config Reader object
         */
        ~ConfigReader();

        /**
         * @brief Create model from configuration file
         *
         * @return Model
         */
        MODEL::Model<T> get_model();

        /**
         * @brief Get number of time steps based on configuration
         *
         * @return Number of time steps
         */
        int get_num_time_steps();

        /**
         * @brief Get the time step size based on configuration
         *
         * @return Time step size
         */
        T get_time_step_size();

        /**
         * @brief Get an integer simulation parameter
         *
         * @param key The key of the simulation parameter
         * @return Value of the simulation parameter
         */
        int get_int_simulation_parameter(std::string key);

        /**
         * @brief Get an integer simulation parameter with default value
         *
         * @param key The key of the simulation parameter
         * @param def Default value if key not in simulation parameters
         * @return Value of the simulation parameter
         */
        int get_int_simulation_parameter(std::string key, int def);

        /**
         * @brief Get a scalar simulation parameter
         *
         * @param key The key of the simulation parameter
         * @return Value of the simulation parameter
         */
        T get_scalar_simulation_parameter(std::string key);

        /**
         * @brief Get a scalar simulation parameter with default value
         *
         * @param key The key of the simulation parameter
         * @param def Default value if key not in simulation parameters
         * @return Value of the simulation parameter
         */
        T get_scalar_simulation_parameter(std::string key, T def);

        /**
         * @brief Get a bool simulation parameter
         *
         * @param key The key of the simulation parameter
         * @return Value of the simulation parameter
         */
        bool get_bool_simulation_parameter(std::string key);

        /**
         * @brief Get a bool simulation parameter with default value
         *
         * @param key The key of the simulation parameter
         * @param def Default value if key not in simulation parameters
         * @return Value of the simulation parameter
         */
        bool get_bool_simulation_parameter(std::string key, bool def);

        T cardiac_cycle_period = 1.0; ///< Cardiac cycle period

    private:
        bool model_created = false;
        simdjson::dom::parser parser;
        simdjson::dom::element config;
        simdjson::dom::element sim_params;
        static bool has_key(simdjson::dom::element &ele, std::string key);
        static T get_default(simdjson::dom::element &ele, std::string key, T def);
        static int get_default(simdjson::dom::element &ele, std::string key, int def);
        static bool get_default(simdjson::dom::element &ele, std::string key, bool def);
        static simdjson::dom::element get_default(simdjson::dom::element &ele, std::string key, simdjson::dom::element def);
        MODEL::TimeDependentParameter<T> get_time_dependent_parameter(simdjson::dom::element &json_times, simdjson::dom::element &json_values);
    };

    template <typename T>
    ConfigReader<T>::ConfigReader(std::string filename)
    {
        config = parser.load(filename);
        sim_params = config["simulation_parameters"];
    }

    template <typename T>
    ConfigReader<T>::~ConfigReader()
    {
    }

    template <typename T>
    bool ConfigReader<T>::has_key(simdjson::dom::element &ele, std::string key)
    {
        // Check if config has key by trying
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
    int ConfigReader<T>::get_default(simdjson::dom::element &ele, std::string key, int def)
    {
        int value = def;
        if (has_key(ele, key))
        {
            value = ele[key].get_int64();
        }
        return value;
    }

    template <typename T>
    bool ConfigReader<T>::get_default(simdjson::dom::element &ele, std::string key, bool def)
    {
        bool value = def;
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
            times.push_back(0.0);
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
        MODEL::TimeDependentParameter<T> param(times, values);
        if (param.isconstant == false)
        {
            cardiac_cycle_period = param.cycle_period;
        }
        return param;
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
                // Create the junction and add to blocks
                std::string junction_name = static_cast<std::string>(junction_config["junction_name"]);
                model.blocks.insert({junction_name, MODEL::Junction<T>(junction_name)});

                // Check for connections to inlet and outlet vessels and append to connections list
                for (simdjson::dom::element inlet_vessel : junction_config["inlet_vessels"])
                {
                    connections.push_back({"V" + std::to_string(inlet_vessel.get_int64()), junction_name});
                }
                for (simdjson::dom::element outlet_vessel : junction_config["outlet_vessels"])
                {
                    connections.push_back({junction_name, "V" + std::to_string(outlet_vessel.get_int64())});
                }
            }
            else
            {
                throw std::invalid_argument("Unknown junction type " + static_cast<std::string>(junction_config["junction_type"]));
            }
        }

        // Create vessels
        std::list<std::string> bc_locations = {"inlet", "outlet"};
        for (simdjson::dom::element vessel_config : config["vessels"])
        {
            simdjson::dom::element vessel_values = vessel_config["zero_d_element_values"];
            std::string vessel_id_string = std::to_string(vessel_config["vessel_id"].get_int64());
            std::string vessel_name = "V" + vessel_id_string;
            if (static_cast<std::string>(vessel_config["zero_d_element_type"]) == "BloodVessel")
            {
                T R = vessel_values["R_poiseuille"], C = get_default(vessel_values, "C", 0.0);
                T L = get_default(vessel_values, "L", 0.0), stenosis_coefficient = get_default(vessel_values, "stenosis_coefficient", 0.0);
                model.blocks.insert({vessel_name, MODEL::BloodVessel<T>(R = R, C = C, L = L, stenosis_coefficient = stenosis_coefficient, vessel_name)});
            }
            else
            {
                throw std::invalid_argument("Unknown vessel type " + static_cast<std::string>(vessel_config["vessel_type"]));
            }

            // Create boundary conditions of blood vessel
            if (has_key(vessel_config, "boundary_conditions") == true)
            {
                simdjson::dom::element vessel_bcs = vessel_config["boundary_conditions"];
                for (auto loc : bc_locations)
                {
                    // Check for BC blocks at inlet or outlet
                    if (has_key(vessel_bcs, loc))
                    {
                        std::string bc_name = "BC" + vessel_id_string + "_" + loc;
                        std::string bc_handle = static_cast<std::string>(vessel_bcs[loc]);
                        for (simdjson::dom::element bc_config : config["boundary_conditions"])
                        {
                            if (static_cast<std::string>(bc_config["bc_name"]) == bc_handle)
                            {
                                simdjson::dom::element bc_values = bc_config["bc_values"];
                                if (static_cast<std::string>(bc_config["bc_type"]) == "RCR")
                                {
                                    T Rp = bc_values["Rp"], C = bc_values["C"], Rd = bc_values["Rd"], Pd = bc_values["Pd"];
                                    model.blocks.insert({bc_name, MODEL::WindkesselBC<T>(Rp = Rp, C = C, Rd = Rd, Pd = Pd, bc_name)});
                                    DEBUG_MSG("Created boundary condition " << bc_name);
                                }
                                else if (static_cast<std::string>(bc_config["bc_type"]) == "FLOW")
                                {
                                    simdjson::dom::element Q_json = bc_values["Q"], t_json = get_default(bc_values, "t", simdjson::dom::element());
                                    auto Q = get_time_dependent_parameter(t_json, Q_json);
                                    model.blocks.insert({bc_name, MODEL::FlowReferenceBC<T>(Q = Q, bc_name)});
                                    DEBUG_MSG("Created boundary condition " << bc_name);
                                }
                                else if (static_cast<std::string>(bc_config["bc_type"]) == "RESISTANCE")
                                {
                                    simdjson::dom::element R_json = bc_values["R"], Pd_json = bc_values["Pd"];
                                    simdjson::dom::element t_json = get_default(bc_values, "t", simdjson::dom::element());
                                    auto R = get_time_dependent_parameter(t_json, R_json);
                                    auto Pd = get_time_dependent_parameter(t_json, Pd_json);
                                    model.blocks.insert({bc_name, MODEL::ResistanceBC<T>(R = R, Pd = Pd, bc_name)});
                                    DEBUG_MSG("Created boundary condition " << bc_name);
                                }
                                else if (static_cast<std::string>(bc_config["bc_type"]) == "PRESSURE")
                                {
                                    simdjson::dom::element P_json = bc_values["P"];
                                    simdjson::dom::element t_json = get_default(bc_values, "t", simdjson::dom::element());
                                    auto P = get_time_dependent_parameter(t_json, P_json);
                                    if (P.isconstant == false)
                                    {
                                        cardiac_cycle_period = P.cycle_period;
                                    }
                                    model.blocks.insert({bc_name, MODEL::PressureReferenceBC<T>(P = P, bc_name)});
                                    DEBUG_MSG("Created boundary condition " << bc_name);
                                }
                                else if (static_cast<std::string>(bc_config["bc_type"]) == "CORONARY")
                                {
                                    T Ra = bc_values["Ra1"], Ram = bc_values["Ra2"], Rv = bc_values["Rv1"], Ca = bc_values["Ca"], Cim = bc_values["Cc"];
                                    simdjson::dom::element Pim_json = bc_values["Pim"], Pv_json = bc_values["P_v"];
                                    simdjson::dom::element t_json = get_default(bc_values, "t", simdjson::dom::element());
                                    auto Pim = get_time_dependent_parameter(t_json, Pim_json);
                                    auto Pv = get_time_dependent_parameter(t_json, Pv_json);
                                    model.blocks.insert({bc_name, MODEL::OpenLoopCoronaryBC<T>(Ra = Ra, Ram = Ram, Rv = Rv, Ca = Ca, Cim = Cim, Pim = Pim, Pv = Pv, bc_name)});
                                    DEBUG_MSG("Created boundary condition " << bc_name);
                                }
                                else
                                {
                                    throw std::invalid_argument("Unknown boundary condition type " + static_cast<std::string>(bc_config["bc_type"]));
                                }
                                break;
                            }
                        }
                        if (loc == "inlet")
                        {
                            connections.push_back({bc_name, vessel_name});
                        }
                        else
                        {
                            connections.push_back({vessel_name, bc_name});
                        }
                        DEBUG_MSG("Found connection " << std::get<0>(connection) << "/" << std::get<1>(connection));
                    }
                }
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

    template <typename T>
    int ConfigReader<T>::get_int_simulation_parameter(std::string key)
    {
        return int(config["simulation_parameters"][key]);
    }

    template <typename T>
    int ConfigReader<T>::get_int_simulation_parameter(std::string key, int def)
    {
        return get_default(sim_params, key, def);
    }

    template <typename T>
    T ConfigReader<T>::get_scalar_simulation_parameter(std::string key)
    {
        return int(config["simulation_parameters"][key]);
    }

    template <typename T>
    T ConfigReader<T>::get_scalar_simulation_parameter(std::string key, T def)
    {
        return get_default(sim_params, key, def);
    }

    template <typename T>
    bool ConfigReader<T>::get_bool_simulation_parameter(std::string key)
    {
        return int(config["simulation_parameters"][key]);
    }

    template <typename T>
    bool ConfigReader<T>::get_bool_simulation_parameter(std::string key, bool def)
    {
        return get_default(sim_params, key, def);
    }

} // namespace IO

#endif // SVZERODSOLVER_IO_CONFIGREADER_HPP_