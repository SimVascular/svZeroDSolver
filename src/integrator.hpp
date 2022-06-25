#ifndef SVZERODSOLVER_INTEGRATOR_H_
#define SVZERODSOLVER_INTEGRATOR_H_

#include <map>
#include <Eigen/Dense>
#include "system.hpp"
#include "model.hpp"

template <typename T>
class State
{
public:
    Eigen::Matrix<T, Eigen::Dynamic, 1> y;
    Eigen::Matrix<T, Eigen::Dynamic, 1> ydot;
    State();
    ~State();
    State(const State &state);

    static State Zero(unsigned int n);
};

template <typename T>
class Integrator
{
private:
    T alpha_m;
    T alpha_f;
    T alpha_m_inv;
    T alpha_f_inv;
    T gamma;
    T gamma_inv;
    T time_step_size;
    T time_step_size_inv;
    T y_dot_coeff;
    T rtol;

    int size;
    Eigen::Matrix<T, Eigen::Dynamic, 1> y_af;
    Eigen::Matrix<T, Eigen::Dynamic, 1> ydot_am;
    Eigen::Matrix<T, Eigen::Dynamic, 1> dy;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> lhs;
    Eigen::Matrix<T, Eigen::Dynamic, 1> res;

    System<T> system;

public:
    Integrator(System<T> &system, T time_step_size, T rho, T atol = 1e-5);
    ~Integrator();
    State<T> step(State<T> &state, T time, Model<T> &model, unsigned int max_iter);
};

template <typename T>
State<T>::State()
{
}

template <typename T>
State<T>::~State()
{
}

template <typename T>
State<T>::State(const State &state)
{
    y = state.y;
    ydot = state.ydot;
}

template <typename T>
State<T> State<T>::Zero(unsigned int n)
{
    static State<T> state;
    state.y = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n);
    state.ydot = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n);
    return state;
}

template <typename T>
Integrator<T>::Integrator(System<T> &system, T time_step_size, T rho, T rtol)
{
    alpha_m = 0.5 * (3.0 - rho) / (1.0 + rho);
    alpha_f = 1.0 / (1.0 + rho);
    alpha_m_inv = alpha_m / 1.0;
    alpha_f_inv = alpha_f / 1.0;
    gamma = 0.5 + alpha_m - alpha_f;
    gamma_inv = 1.0 / gamma;

    this->system = system;
    this->time_step_size = time_step_size;
    this->rtol = rtol;
    time_step_size_inv = 1.0 / time_step_size;

    size = system.F.row(0).size();
    y_af = Eigen::Matrix<T, Eigen::Dynamic, 1>(size);
    ydot_am = Eigen::Matrix<T, Eigen::Dynamic, 1>(size);
    dy = Eigen::Matrix<T, Eigen::Dynamic, 1>(size);
    lhs = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(size, size);
    res = Eigen::Matrix<T, Eigen::Dynamic, 1>(size);

    y_dot_coeff = alpha_m / (alpha_f * gamma) * time_step_size_inv;
}

template <typename T>
Integrator<T>::~Integrator()
{
}

template <typename T>
State<T> Integrator<T>::step(State<T> &state, T time, Model<T> &model, unsigned int max_iter)
{
    State<T> old_state = state;
    State<T> new_state;

    new_state.y = old_state.y + 0.5 * time_step_size * old_state.ydot;
    new_state.ydot = old_state.ydot * (gamma - 0.5) * gamma_inv;

    y_af = old_state.y + alpha_f * (new_state.y - old_state.y);
    ydot_am = old_state.ydot + alpha_m * (new_state.ydot - old_state.ydot);

    T new_time = time + alpha_f * time_step_size;

    // Update time in blocks
    model.update_time(system, new_time);

    dy.setZero();
    for (size_t i = 0; i < max_iter; i++)
    {
        // Update solution and assemble
        model.update_solution(system, y_af);

        // Calculate RHS and LHS
        res = -(system.E * ydot_am) - (system.F * y_af) - system.C;
        if (res.cwiseAbs().maxCoeff() < rtol)
        {
            break;
        }
        else if (i == max_iter - 1)
        {
            throw std::runtime_error("Maxium number of non-linear iterations reached.");
        }
        lhs = system.F + system.dE + system.dF + system.dC + system.E * y_dot_coeff;

        // Solve system
        // TODO: Works only if matrix is invertable: Check if True otherwise use colPivHouseholderQr
        dy = lhs.partialPivLu().solve(res);

        // Update solution
        y_af += dy;
        ydot_am += dy * y_dot_coeff;
    }

    new_state.y = old_state.y + (y_af - old_state.y) * alpha_f_inv;
    new_state.ydot = old_state.ydot + (ydot_am - old_state.ydot) / alpha_m_inv;

    return new_state;
}

#endif // SVZERODSOLVER_INTEGRATOR_H_