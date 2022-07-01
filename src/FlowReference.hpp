#ifndef SVZERODSOLVER_FLOWREFERENCE_H_
#define SVZERODSOLVER_FLOWREFERENCE_H_

#include "block.hpp"
#include "system.hpp"
#include "parameter.hpp"

template <typename T>
class FlowReference : public Block<T>
{
public:
    struct Parameters : public Block<T>::Parameters
    {
        TimeDependentParameter<T> Q; // Flow at timestep
    };
    FlowReference(TimeDependentParameter<T> Q, std::string name);
    ~FlowReference();
    void setup_dofs(DOFHandler &dofhandler);

    // Dense
    void update_constant(System<T> &system);
    void update_time(System<T> &system, T time);

    // Sparse
    void update_constant(SparseSystem<T> &system);
    void update_time(SparseSystem<T> &system, T time);

    // Number of triplets that will be added to global matrices (relevant for sparse reservation)
    std::map<std::string, int> num_triplets = {
        {"F", 1},
        {"E", 0},
        {"D", 0},
    };

    void to_steady();

private:
    Parameters params;
};

template <typename T>
FlowReference<T>::FlowReference(TimeDependentParameter<T> Q, std::string name) : Block<T>(name)
{
    this->name = name;
    this->params.Q = Q;
}

template <typename T>
FlowReference<T>::~FlowReference()
{
}

template <typename T>
void FlowReference<T>::setup_dofs(DOFHandler &dofhandler)
{
    Block<T>::setup_dofs_(dofhandler, 1, 0);
}

template <typename T>
void FlowReference<T>::update_constant(System<T> &system)
{
    system.F(this->global_eqn_ids[0], this->global_var_ids[1]) = 1.0;
}

template <typename T>
void FlowReference<T>::update_time(System<T> &system, T time)
{
    system.C(this->global_eqn_ids[0]) = -params.Q.get(time);
}

template <typename T>
void FlowReference<T>::update_constant(SparseSystem<T> &system)
{
    system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) = 1.0;
}

template <typename T>
void FlowReference<T>::update_time(SparseSystem<T> &system, T time)
{
    system.C(this->global_eqn_ids[0]) = -params.Q.get(time);
}

template <typename T>
void FlowReference<T>::to_steady()
{
    params.Q.to_steady();
}

#endif // SVZERODSOLVER_FLOWREFERENCE_H_