#ifndef SVZERODSOLVER_MODEL_H_
#define SVZERODSOLVER_MODEL_H_

#include <map>
#include <list>

#include "junction.hpp"
#include "bloodvessel.hpp"
#include "rcrblockwithdistalpressure.hpp"
#include "flowreference.hpp"
#include "dofhandler.hpp"
#include "node.hpp"
#include "system.hpp"
#include "parameter.hpp"

template <typename T>
class Model
{
public:
    Model();
    ~Model();

    std::map<std::string, std::variant<Junction<T>, BloodVessel<T>, FlowReference<T>, RCRBlockWithDistalPressure<T>>> blocks;
    DOFHandler dofhandler;
    std::list<Node> nodes;
    std::list<TimeDependentParameter<T>> time_params;

    // Dense
    void update_constant(System<T> &system);
    void update_time(System<T> &system, T time);
    void update_solution(System<T> &system, Eigen::Matrix<T, Eigen::Dynamic, 1> &y);

    // Sparse
    void update_constant(SparseSystem<T> &system);
    void update_time(SparseSystem<T> &system, T time);
    void update_solution(SparseSystem<T> &system, Eigen::Matrix<T, Eigen::Dynamic, 1> &y);

    void to_steady();
    std::map<std::string, int> get_num_triplets();
};

template <typename T>
Model<T>::Model()
{
}

template <typename T>
Model<T>::~Model()
{
}

template <typename T>
void Model<T>::update_constant(System<T> &system)
{
    for (auto &&elem : blocks)
    {
        std::visit([&](auto &&block)
                   { block.update_constant(system); },
                   elem.second);
    }
}

template <typename T>
void Model<T>::update_time(System<T> &system, T time)
{
    for (auto &&elem : blocks)
    {
        std::visit([&](auto &&block)
                   { block.update_time(system, time); },
                   elem.second);
    }
}

template <typename T>
void Model<T>::update_solution(System<T> &system, Eigen::Matrix<T, Eigen::Dynamic, 1> &y)
{
    for (auto &&elem : blocks)
    {
        std::visit([&](auto &&block)
                   { block.update_solution(system, y); },
                   elem.second);
    }
}

template <typename T>
void Model<T>::update_constant(SparseSystem<T> &system)
{
    for (auto &&elem : blocks)
    {
        std::visit([&](auto &&block)
                   { block.update_constant(system); },
                   elem.second);
    }
}

template <typename T>
void Model<T>::update_time(SparseSystem<T> &system, T time)
{
    for (auto &&elem : blocks)
    {
        std::visit([&](auto &&block)
                   { block.update_time(system, time); },
                   elem.second);
    }
}

template <typename T>
void Model<T>::update_solution(SparseSystem<T> &system, Eigen::Matrix<T, Eigen::Dynamic, 1> &y)
{
    for (auto &&elem : blocks)
    {
        std::visit([&](auto &&block)
                   { block.update_solution(system, y); },
                   elem.second);
    }
}

template <typename T>
void Model<T>::to_steady()
{
    for (auto &&elem : blocks)
    {
        std::visit([&](auto &&block)
                   { block.to_steady(); },
                   elem.second);
    }
}

template <typename T>
std::map<std::string, int> Model<T>::get_num_triplets()
{
    std::map<std::string, int> num_triplets = {
        {"F", 0},
        {"E", 0},
        {"dF", 0},
        {"dE", 0},
        {"dC", 0},
    };
    for (auto &&elem : blocks)
    {
        std::visit([&](auto &&block)
                   { for (auto &[key, value] : block.num_triplets){num_triplets[key] += value;} },
                   elem.second);
    }
    return num_triplets;
}

#endif // SVZERODSOLVER_MODEL_H_