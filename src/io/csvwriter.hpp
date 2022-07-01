#ifndef SVZERODSOLVER_IO_CSVWRITER_HPP_
#define SVZERODSOLVER_IO_CSVWRITER_HPP_

#include <string>
#include <vector>
#include <fstream>

#include "../algebra/state.hpp"
#include "../model/model.hpp"
#include "../helpers/startswith.hpp"

namespace IO
{

    template <typename T>
    void write_csv(std::string path, std::vector<T> times, std::vector<ALGEBRA::State<T>> states, MODEL::Model<T> model, bool mean = false)
    {

        std::stringstream out;
        out << "name,time,flow_in,flow_out,pressure_in,pressure_out\n";
        char buff[100];
        T num_steps = times.size();

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
                if (mean)
                {
                    T inflow_mean = 0.0;
                    T outflow_mean = 0.0;
                    T inpres_mean = 0.0;
                    T outpres_mean = 0.0;
                    for (size_t i = 0; i < times.size(); i++)
                    {
                        inflow_mean += states[i].y[inflow_dof];
                        outflow_mean += states[i].y[outflow_dof];
                        inpres_mean += states[i].y[inpres_dof];
                        outpres_mean += states[i].y[outpres_dof];
                    }
                    inflow_mean /= num_steps;
                    outflow_mean /= num_steps;
                    inpres_mean /= num_steps;
                    outpres_mean /= num_steps;
                    sprintf(buff, "%s,,%.10e,%.10e,%.10e,%.10e\n", name.c_str(), inflow_mean, outflow_mean, inpres_mean, outpres_mean);
                    out << buff;
                }
                else
                {
                    for (size_t i = 0; i < times.size(); i++)
                    {
                        sprintf(buff, "%s,%.10f,%.10e,%.10e,%.10e,%.10e\n", name.c_str(), times[i], states[i].y[inflow_dof], states[i].y[outflow_dof], states[i].y[inpres_dof], states[i].y[outpres_dof]);
                        out << buff;
                    }
                }
            }
        }

        std::ofstream ofs(path);
        ofs << out.str();
        ofs.close();
    }
} // namespace IO

#endif // SVZERODSOLVER_IO_CSVWRITER_HPP_