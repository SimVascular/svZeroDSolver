#ifndef SVZERODSOLVER_BLOODVESSEL_H_
#define SVZERODSOLVER_BLOODVESSEL_H_

#include "block.hpp"
#include <math.h>

template <typename T>
class BloodVessel : public Block<T>
{
public:
    struct Parameters : public Block<T>::Parameters
    {
        T R;                    // Poseuille resistance
        T C;                    // Capacitance
        T L;                    // Inductance
        T stenosis_coefficient; // Stenosis Coefficient
    };
    BloodVessel(T R, T C, T L, T stenosis_coefficient, std::string name);
    ~BloodVessel();
    void setup_dofs(DOFHandler &dofhandler);

    // Dense
    void update_constant(System<T> &system);
    void update_solution(System<T> &system, Eigen::Matrix<T, Eigen::Dynamic, 1> &y);

    // Sparse
    void update_constant(SparseSystem<T> &system);
    void update_solution(SparseSystem<T> &system, Eigen::Matrix<T, Eigen::Dynamic, 1> &y);

    // Number of triplets that will be added to global matrices (relevant for sparse reservation)
    std::map<std::string, int> num_triplets = {
        {"F", 10},
        {"E", 2},
        {"dF", 2},
        {"dE", 0},
        {"dC", 0},
    };

    std::string name;

private:
    Parameters params;
};

template <typename T>
BloodVessel<T>::BloodVessel(T R, T C, T L, T stenosis_coefficient, std::string name) : Block<T>(name)
{
    this->name = name;
    this->params.R = R;
    this->params.C = C;
    this->params.L = L;
    this->params.stenosis_coefficient = stenosis_coefficient;
}

template <typename T>
BloodVessel<T>::~BloodVessel()
{
}

template <typename T>
void BloodVessel<T>::setup_dofs(DOFHandler &dofhandler)
{
    Block<T>::setup_dofs_(dofhandler, 3, 1);
}

template <typename T>
void BloodVessel<T>::update_constant(System<T> &system)
{
    system.E(this->global_eqn_ids[0], this->global_var_ids[3]) = -params.L;
    system.E(this->global_eqn_ids[1], this->global_var_ids[4]) = -params.C;

    system.F(this->global_eqn_ids[0], this->global_var_ids[0]) = 1.0;
    system.F(this->global_eqn_ids[0], this->global_var_ids[1]) = -params.R;
    system.F(this->global_eqn_ids[0], this->global_var_ids[2]) = -1.0;

    system.F(this->global_eqn_ids[1], this->global_var_ids[1]) = 1.0;
    system.F(this->global_eqn_ids[1], this->global_var_ids[3]) = -1.0;

    system.F(this->global_eqn_ids[2], this->global_var_ids[0]) = 1.0;
    system.F(this->global_eqn_ids[2], this->global_var_ids[1]) = -params.R;
    system.F(this->global_eqn_ids[2], this->global_var_ids[4]) = -1.0;
}

template <typename T>
void BloodVessel<T>::update_solution(System<T> &system, Eigen::Matrix<T, Eigen::Dynamic, 1> &y)
{
    T q_in = fabs(y[this->inlet_nodes[0]->flow_dof]);
    T fac1 = -params.stenosis_coefficient * q_in;
    T fac2 = fac1 - params.R;
    system.F(this->global_eqn_ids[0], this->global_var_ids[1]) = fac2;
    system.F(this->global_eqn_ids[2], this->global_var_ids[1]) = fac2;
    system.dF(this->global_eqn_ids[0], this->global_var_ids[1]) = fac1;
    system.dF(this->global_eqn_ids[2], this->global_var_ids[1]) = fac1;
}

template <typename T>
void BloodVessel<T>::update_constant(SparseSystem<T> &system)
{
    system.E.coeffRef(this->global_eqn_ids[0], this->global_var_ids[3]) = -params.L;
    system.E.coeffRef(this->global_eqn_ids[1], this->global_var_ids[4]) = -params.C;

    system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[0]) = 1.0;
    system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) = -params.R;
    system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[2]) = -1.0;

    system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[1]) = 1.0;
    system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[3]) = -1.0;

    system.F.coeffRef(this->global_eqn_ids[2], this->global_var_ids[0]) = 1.0;
    system.F.coeffRef(this->global_eqn_ids[2], this->global_var_ids[1]) = -params.R;
    system.F.coeffRef(this->global_eqn_ids[2], this->global_var_ids[4]) = -1.0;
}

template <typename T>
void BloodVessel<T>::update_solution(SparseSystem<T> &system, Eigen::Matrix<T, Eigen::Dynamic, 1> &y)
{
    T q_in = fabs(y[this->inlet_nodes[0]->flow_dof]);
    T fac1 = -params.stenosis_coefficient * q_in;
    T fac2 = fac1 - params.R;
    system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) = fac2;
    system.F.coeffRef(this->global_eqn_ids[2], this->global_var_ids[1]) = fac2;
    system.dF.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) = fac1;
    system.dF.coeffRef(this->global_eqn_ids[2], this->global_var_ids[1]) = fac1;
}

#endif // SVZERODSOLVER_BLOODVESSEL_H_