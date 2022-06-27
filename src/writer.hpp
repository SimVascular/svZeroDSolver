#ifndef SVZERODSOLVER_WRITER_H_
#define SVZERODSOLVER_WRITER_H_

#include "json.h"
#include <string>
#include <vector>
#include "integrator.hpp"
#include "model.hpp"

bool startsWith(const std::string &str, const std::string &prefix)
{
    return str.size() >= prefix.size() && str.compare(0, prefix.size(), prefix) == 0;
}

template <typename T>
void write_json(std::string path, std::vector<T> times, std::vector<State<T>> states, Model<T> model)
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
                   { if (startsWith(block.name, "V")){name = block.name; inflow_dof = block.inlet_nodes[0]->flow_dof; outflow_dof = block.outlet_nodes[0]->flow_dof; inpres_dof = block.inlet_nodes[0]->pres_dof; outpres_dof = block.outlet_nodes[0]->pres_dof;} },
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

template <typename T>
void write_csv(std::string path, std::vector<T> times, std::vector<State<T>> states, Model<T> model)
{

    std::stringstream out;
    out << "name,time,flow_in,flow_out,pressure_in,pressure_out\n";
    char buff[100];

    for (auto &[key, elem] : model.blocks)
    {
        std::string name = "NoName";
        unsigned int inflow_dof;
        unsigned int outflow_dof;
        unsigned int inpres_dof;
        unsigned int outpres_dof;
        std::visit([&](auto &&block)
                   { if (startsWith(block.name, "V")){name = block.name; inflow_dof = block.inlet_nodes[0]->flow_dof; outflow_dof = block.outlet_nodes[0]->flow_dof; inpres_dof = block.inlet_nodes[0]->pres_dof; outpres_dof = block.outlet_nodes[0]->pres_dof;} },
                   elem);

        if (name != "NoName")
        {
            for (size_t i = 0; i < times.size(); i++)
            {
                // out << name << "," << times[i] << "," << states[i].y[inflow_dof] << "," << states[i].y[outflow_dof] << "," << states[i].y[inpres_dof] << "," << states[i].y[outpres_dof] << "\n";
                sprintf(buff, "%s,%.10f,%.10e,%.10e,%.10e,%.10e\n", name.c_str(), times[i], states[i].y[inflow_dof], states[i].y[outflow_dof], states[i].y[inpres_dof], states[i].y[outpres_dof]);
                out << buff;
            }
        }
    }

    std::ofstream ofs(path);
    ofs << out.str();
    ofs.close();
}

#endif // SVZERODSOLVER_WRITER_H_