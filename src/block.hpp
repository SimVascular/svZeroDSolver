#ifndef SVZERODSOLVER_BLOCK_H_
#define SVZERODSOLVER_BLOCK_H_

#include <vector>
#include <map>
#include "node.hpp"
#include "dofhandler.hpp"
#include "system.hpp"

template <typename T>
class Block
{
public:
    struct Parameters
    {
    };
    Block();
    Block(std::string name);
    ~Block();

    std::string name;
    std::vector<Node *> inlet_nodes;
    std::vector<Node *> outlet_nodes;

    std::vector<unsigned int> global_var_ids;
    std::vector<unsigned int> global_eqn_ids;

    void setup_dofs_(DOFHandler &dofhandler, unsigned int num_equations, unsigned int num_internal_vars);

    // Dense
    void update_constant(System<T> &system);
    void update_time(System<T> &system, T time);
    void update_solution(System<T> &system, Eigen::Matrix<T, Eigen::Dynamic, 1> &y);

    // Sparse
    void update_constant(SparseSystem<T> &system);
    void update_time(SparseSystem<T> &system, T time);
    void update_solution(SparseSystem<T> &system, Eigen::Matrix<T, Eigen::Dynamic, 1> &y);

    // Number of triplets that will be added to global matrices (relevant for sparse reservation)
    std::map<std::string, int> num_triplets = {
        {"F", 0},
        {"E", 0},
        {"D", 0},
    };

    void to_steady();

private:
    Parameters params;
};

template <typename T>
Block<T>::Block()
{
}

template <typename T>
Block<T>::Block(std::string name)
{
    this->name = name;
}

template <typename T>
Block<T>::~Block()
{
}

template <typename T>
void Block<T>::setup_dofs_(DOFHandler &dofhandler, unsigned int num_equations, unsigned int num_internal_vars)
{
    // Collect external DOFs from inlet and outlet nodes
    for (auto inlet_node : inlet_nodes)
    {
        global_var_ids.push_back(inlet_node->pres_dof);
        global_var_ids.push_back(inlet_node->flow_dof);
    }
    for (auto outlet_node : outlet_nodes)
    {
        global_var_ids.push_back(outlet_node->pres_dof);
        global_var_ids.push_back(outlet_node->flow_dof);
    }

    // Register internal variables of block
    for (unsigned int i = 0; i < num_internal_vars; i++)
    {
        global_var_ids.push_back(dofhandler.register_variable("var_" + std::to_string(i) + "_" + name));
    }

    // Register equations of block
    for (unsigned int i = 0; i < num_equations; i++)
    {
        global_eqn_ids.push_back(dofhandler.register_equation());
    }
}

template <typename T>
void Block<T>::update_constant(System<T> &system)
{
}

template <typename T>
void Block<T>::update_time(System<T> &system, T time)
{
}

template <typename T>
void Block<T>::update_solution(System<T> &system, Eigen::Matrix<T, Eigen::Dynamic, 1> &y)
{
}

template <typename T>
void Block<T>::update_constant(SparseSystem<T> &system)
{
}

template <typename T>
void Block<T>::update_time(SparseSystem<T> &system, T time)
{
}

template <typename T>
void Block<T>::update_solution(SparseSystem<T> &system, Eigen::Matrix<T, Eigen::Dynamic, 1> &y)
{
}

template <typename T>
void Block<T>::to_steady()
{
}

#endif // SVZERODSOLVER_BLOCK_H_