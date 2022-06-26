#ifndef SVZERODSOLVER_JUNCTION_H_
#define SVZERODSOLVER_JUNCTION_H_

#include "block.hpp"

template <typename T>
class Junction : public Block<T>
{
public:
    struct Parameters : public Block<T>::Parameters
    {
    };
    Junction(std::string name);
    ~Junction();
    void setup_dofs(DOFHandler &dofhandler);

    // Dense
    void update_constant(System<T> &system);

    // Sparse
    void update_constant(SparseSystem<T> &system);

private:
    Parameters params;
    unsigned int num_inlets;
    unsigned int num_outlets;
};

template <typename T>
Junction<T>::Junction(std::string name) : Block<T>(name)
{
    this->name = name;
}

template <typename T>
Junction<T>::~Junction()
{
}

template <typename T>
void Junction<T>::setup_dofs(DOFHandler &dofhandler)
{
    // Set number of equations of a junction block based on number of
    // inlets/outlets. Must be set before calling parent constructor
    num_inlets = this->inlet_nodes.size();
    num_outlets = this->outlet_nodes.size();
    Block<T>::setup_dofs_(dofhandler, num_inlets + num_outlets, 0);
}

template <typename T>
void Junction<T>::update_constant(System<T> &system)
{
    for (size_t i = 0; i < (num_inlets + num_outlets - 1); i++)
    {
        system.F(this->global_eqn_ids[i], this->global_var_ids[0]) = 1.0;
        system.F(this->global_eqn_ids[i], this->global_var_ids[2 * i + 2]) = -1.0;
    }
    for (size_t i = 1; i < num_inlets * 2; i = i + 2)
    {
        system.F(this->global_eqn_ids[num_inlets + num_outlets - 1], this->global_var_ids[i]) = 1.0;
    }
    for (size_t i = (num_inlets * 2) + 1; i < (num_inlets + num_outlets) * 2; i = i + 2)
    {
        system.F(this->global_eqn_ids[num_inlets + num_outlets - 1], this->global_var_ids[i]) = -1.0;
    }
}

template <typename T>
void Junction<T>::update_constant(SparseSystem<T> &system)
{
    for (size_t i = 0; i < (num_inlets + num_outlets - 1); i++)
    {
        system.F.coeffRef(this->global_eqn_ids[i], this->global_var_ids[0]) = 1.0;
        system.F.coeffRef(this->global_eqn_ids[i], this->global_var_ids[2 * i + 2]) = -1.0;
    }
    for (size_t i = 1; i < num_inlets * 2; i = i + 2)
    {
        system.F.coeffRef(this->global_eqn_ids[num_inlets + num_outlets - 1], this->global_var_ids[i]) = 1.0;
    }
    for (size_t i = (num_inlets * 2) + 1; i < (num_inlets + num_outlets) * 2; i = i + 2)
    {
        system.F.coeffRef(this->global_eqn_ids[num_inlets + num_outlets - 1], this->global_var_ids[i]) = -1.0;
    }
}

#endif // SVZERODSOLVER_JUNCTION_H_