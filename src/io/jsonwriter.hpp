/**
 * @file jsonwriter.hpp
 * @brief IO::write_json source file
 */
#ifndef SVZERODSOLVER_IO_JSONWRITER_HPP_
#define SVZERODSOLVER_IO_JSONWRITER_HPP_

#include <json.h>
#include <vector>
#include "../algebra/state.hpp"
#include "../model/model.hpp"
#include "../helpers/startswith.hpp"

namespace IO
{

    /**
     * @brief Write the solution to a json file
     *
     * @tparam T Scalar type (e.g. `float`, `double`)
     * @param path Path to the output json file
     * @param times Sequence of time steps corresponding to the solutions
     * @param states Sequence of states corresponding to the time steps
     * @param model The underlying model
     */
    template <typename T>
    void write_json(std::string path, std::vector<T> times, std::vector<ALGEBRA::State<T>> states, MODEL::Model<T> model)
    {
        Json::Value output;
        Json::Value json_times(Json::arrayValue);
        for (auto time : times)
        {
            json_times.append(Json::Value(time));
        }

        Json::Value json_names(Json::arrayValue);
        Json::Value json_flow_in(Json::arrayValue);
        Json::Value json_flow_out(Json::arrayValue);
        Json::Value json_pres_in(Json::arrayValue);
        Json::Value json_pres_out(Json::arrayValue);

        for (auto &[key, elem] : model.blocks)
        {
            std::string name = "NoName";
            unsigned int inflow_dof;
            unsigned int outflow_dof;
            unsigned int inpres_dof;
            unsigned int outpres_dof;
            std::visit([&](auto &&block)
                       { if (HELPERS::startswith(block.name, "V")){name = block.name; inflow_dof = block.inlet_nodes[0]->flow_dof; outflow_dof = block.outlet_nodes[0]->flow_dof; inpres_dof = block.inlet_nodes[0]->pres_dof; outpres_dof = block.outlet_nodes[0]->pres_dof;} },
                       elem);

            if (name != "NoName")
            {
                json_names.append(name);
                Json::Value json_flow_in_i(Json::arrayValue);
                Json::Value json_flow_out_i(Json::arrayValue);
                Json::Value json_pres_in_i(Json::arrayValue);
                Json::Value json_pres_out_i(Json::arrayValue);
                for (auto state : states)
                {
                    json_flow_in_i.append(state.y[inflow_dof]);
                    json_flow_out_i.append(state.y[outflow_dof]);
                    json_pres_in_i.append(state.y[inpres_dof]);
                    json_pres_out_i.append(state.y[outpres_dof]);
                }
                json_flow_in.append(json_flow_in_i);
                json_flow_out.append(json_flow_out_i);
                json_pres_in.append(json_pres_in_i);
                json_pres_out.append(json_pres_out_i);
            }
        }

        output["time"] = json_times;
        output["names"] = json_names;
        output["flow_in"] = json_flow_in;
        output["flow_out"] = json_flow_out;
        output["pressure_in"] = json_pres_in;
        output["pressure_out"] = json_pres_out;

        Json::FastWriter writer;
        std::ofstream out(path);
        out << writer.write(output);
        out.close();
    }
} // namespace IO

#endif // SVZERODSOLVER_IO_JSONWRITER_HPP_