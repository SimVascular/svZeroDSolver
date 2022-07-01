#ifndef SVZERODSOLVER_MODEL_MODEL_HPP_
#define SVZERODSOLVER_MODEL_MODEL_HPP_

#include <map>
#include <list>

#include "../algebra/densesystem.hpp"
#include "../algebra/sparsesystem.hpp"

#include "junction.hpp"
#include "bloodvessel.hpp"
#include "rcrblockwithdistalpressure.hpp"
#include "flowreference.hpp"
#include "dofhandler.hpp"
#include "node.hpp"

namespace MODEL
{

    template <typename T>
    class Model
    {
    public:
        Model();
        ~Model();

        std::map<std::string, std::variant<Junction<T>, BloodVessel<T>, FlowReference<T>, RCRBlockWithDistalPressure<T>>> blocks;
        DOFHandler dofhandler;
        std::list<Node> nodes;

        // Dense
        void update_constant(ALGEBRA::DenseSystem<T> &system);
        void update_time(ALGEBRA::DenseSystem<T> &system, T time);
        void update_solution(ALGEBRA::DenseSystem<T> &system, Eigen::Matrix<T, Eigen::Dynamic, 1> &y);

        // Sparse
        void update_constant(ALGEBRA::SparseSystem<T> &system);
        void update_time(ALGEBRA::SparseSystem<T> &system, T time);
        void update_solution(ALGEBRA::SparseSystem<T> &system, Eigen::Matrix<T, Eigen::Dynamic, 1> &y);

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
    void Model<T>::update_constant(ALGEBRA::DenseSystem<T> &system)
    {
        for (auto &&elem : blocks)
        {
            std::visit([&](auto &&block)
                       { block.update_constant(system); },
                       elem.second);
        }
    }

    template <typename T>
    void Model<T>::update_time(ALGEBRA::DenseSystem<T> &system, T time)
    {
        for (auto &&elem : blocks)
        {
            std::visit([&](auto &&block)
                       { block.update_time(system, time); },
                       elem.second);
        }
    }

    template <typename T>
    void Model<T>::update_solution(ALGEBRA::DenseSystem<T> &system, Eigen::Matrix<T, Eigen::Dynamic, 1> &y)
    {
        for (auto &&elem : blocks)
        {
            std::visit([&](auto &&block)
                       { block.update_solution(system, y); },
                       elem.second);
        }
    }

    template <typename T>
    void Model<T>::update_constant(ALGEBRA::SparseSystem<T> &system)
    {
        for (auto &&elem : blocks)
        {
            std::visit([&](auto &&block)
                       { block.update_constant(system); },
                       elem.second);
        }
    }

    template <typename T>
    void Model<T>::update_time(ALGEBRA::SparseSystem<T> &system, T time)
    {
        for (auto &&elem : blocks)
        {
            std::visit([&](auto &&block)
                       { block.update_time(system, time); },
                       elem.second);
        }
    }

    template <typename T>
    void Model<T>::update_solution(ALGEBRA::SparseSystem<T> &system, Eigen::Matrix<T, Eigen::Dynamic, 1> &y)
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
            {"D", 0},
        };
        for (auto &&elem : blocks)
        {
            std::visit([&](auto &&block)
                       { for (auto &[key, value] : block.num_triplets){num_triplets[key] += value;} },
                       elem.second);
        }
        return num_triplets;
    }

} // namespace MODEL

#endif // SVZERODSOLVER_MODEL_MODEL_HPP_