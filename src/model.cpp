#include "model.hpp"

Model::Model()
{
}
Model::~Model()
{
}

void Model::update_constant(System &system)
{
    for (auto &&elem : blocks)
    {
        std::visit([&](auto &&block)
                   { block.update_constant(system); },
                   elem.second);
    }
}
void Model::update_time(System &system, double time)
{
    for (auto &&elem : blocks)
    {
        std::visit([&](auto &&block)
                   { block.update_time(system, time); },
                   elem.second);
    }
}
void Model::update_solution(System &system, Eigen::VectorXd &y)
{
    for (auto &&elem : blocks)
    {
        std::visit([&](auto &&block)
                   { block.update_solution(system, y); },
                   elem.second);
    }
}