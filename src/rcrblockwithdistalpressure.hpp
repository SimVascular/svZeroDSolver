#ifndef SVZERODSOLVER_RCRBLOCKWITHDISTALPRESSURE_H_
#define SVZERODSOLVER_RCRBLOCKWITHDISTALPRESSURE_H_

#include "block.hpp"

template <typename T>
class RCRBlockWithDistalPressure : public Block<T>
{
public:
    struct Parameters : public Block<T>::Parameters
    {
        T Rp; // Proximal resistance
        T C;  // Capacitance
        T Rd; // Distal restistance
        T Pd; // Distal Pressure
    };
    RCRBlockWithDistalPressure(T Rp, T C, T Rd, T Pd, std::string name);
    ~RCRBlockWithDistalPressure();
    void setup_dofs(DOFHandler &dofhandler);

    // Dense
    void update_constant(System<T> &system);
    void update_time(System<T> &system, T time);

    // Sparse
    void update_constant(SparseSystem<T> &system);
    void update_time(SparseSystem<T> &system, T time);

    // Number of triplets that will be added to global matrices (relevant for sparse reservation)
    std::map<std::string, int> num_triplets = {
        {"F", 5},
        {"E", 1},
        {"dF", 0},
        {"dE", 0},
        {"dC", 0},
    };

    void to_steady();

private:
    Parameters params;
};

template <typename T>
RCRBlockWithDistalPressure<T>::RCRBlockWithDistalPressure(T Rp, T C, T Rd, T Pd, std::string name) : Block<T>(name)
{
    this->name = name;
    this->params.Rp = Rp;
    this->params.C = C;
    this->params.Rd = Rd;
    this->params.Pd = Pd;
}

template <typename T>
RCRBlockWithDistalPressure<T>::~RCRBlockWithDistalPressure()
{
}

template <typename T>
void RCRBlockWithDistalPressure<T>::setup_dofs(DOFHandler &dofhandler)
{
    Block<T>::setup_dofs_(dofhandler, 2, 1);
}

template <typename T>
void RCRBlockWithDistalPressure<T>::update_constant(System<T> &system)
{
    system.F(this->global_eqn_ids[0], this->global_var_ids[0]) = 1.0;
    system.F(this->global_eqn_ids[0], this->global_var_ids[2]) = -1.0;
    system.F(this->global_eqn_ids[1], this->global_var_ids[2]) = -1.0;
}

template <typename T>
void RCRBlockWithDistalPressure<T>::update_time(System<T> &system, T time)
{
    system.E(this->global_eqn_ids[1], this->global_var_ids[2]) = -params.Rd * params.C;
    system.F(this->global_eqn_ids[0], this->global_var_ids[1]) = -params.Rp;
    system.F(this->global_eqn_ids[1], this->global_var_ids[1]) = params.Rd;
    system.C(this->global_eqn_ids[1]) = params.Pd;
}

template <typename T>
void RCRBlockWithDistalPressure<T>::update_constant(SparseSystem<T> &system)
{
    system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[0]) = 1.0;
    system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[2]) = -1.0;
    system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[2]) = -1.0;
}

template <typename T>
void RCRBlockWithDistalPressure<T>::update_time(SparseSystem<T> &system, T time)
{
    system.E.coeffRef(this->global_eqn_ids[1], this->global_var_ids[2]) = -params.Rd * params.C;
    system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) = -params.Rp;
    system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[1]) = params.Rd;
    system.C(this->global_eqn_ids[1]) = params.Pd;
}

template <typename T>
void RCRBlockWithDistalPressure<T>::to_steady()
{
    params.C = 0.0;
}

#endif // SVZERODSOLVER_RCRBLOCKWITHDISTALPRESSURE_H_