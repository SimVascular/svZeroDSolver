#ifndef SVZERODSOLVER_ALGEBRA_INTEGRATOR_HPP_
#define SVZERODSOLVER_ALGEBRA_INTEGRATOR_HPP_

#include <map>

#include "../model/model.hpp"
#include "../external/eigen/Eigen/Dense"
#include "state.hpp"

namespace ALGEBRA
{

    template <typename T, template <class> class S>
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
        T atol;

        int size;
        Eigen::Matrix<T, Eigen::Dynamic, 1> y_af;
        Eigen::Matrix<T, Eigen::Dynamic, 1> ydot_am;

        S<T> system;

    public:
        Integrator(S<T> &system, T time_step_size, T rho, T atol = 1e-5);
        ~Integrator();
        State<T> step(State<T> &state, T time, MODEL::Model<T> &model, unsigned int max_iter);
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

    template <typename T, template <class> class S>
    Integrator<T, S>::Integrator(S<T> &system, T time_step_size, T rho, T atol)
    {
        alpha_m = 0.5 * (3.0 - rho) / (1.0 + rho);
        alpha_f = 1.0 / (1.0 + rho);
        alpha_m_inv = alpha_m / 1.0;
        alpha_f_inv = alpha_f / 1.0;
        gamma = 0.5 + alpha_m - alpha_f;
        gamma_inv = 1.0 / gamma;

        this->system = system;
        this->time_step_size = time_step_size;
        this->atol = atol;
        time_step_size_inv = 1.0 / time_step_size;

        size = system.F.row(0).size();
        y_af = Eigen::Matrix<T, Eigen::Dynamic, 1>(size);
        ydot_am = Eigen::Matrix<T, Eigen::Dynamic, 1>(size);

        y_dot_coeff = alpha_m / (alpha_f * gamma) * time_step_size_inv;
    }

    template <typename T, template <class> class S>
    Integrator<T, S>::~Integrator()
    {
    }

    template <typename T, template <class> class S>
    State<T> Integrator<T, S>::step(State<T> &state, T time, MODEL::Model<T> &model, unsigned int max_iter)
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

        for (size_t i = 0; i < max_iter; i++)
        {
            // Update solution and assemble
            model.update_solution(system, y_af);

            // Determine residuum and check termination criteria
            system.update_residual(y_af, ydot_am);
            if (system.residual.cwiseAbs().maxCoeff() < atol)
            {
                break;
            }
            else if (i == max_iter - 1)
            {
                throw std::runtime_error("Maxium number of non-linear iterations reached.");
            }

            // Determine jacobian
            system.update_jacobian(y_dot_coeff);

            // Solve system
            system.solve();

            // Update solution
            y_af += system.dy;
            ydot_am += system.dy * y_dot_coeff;
        }

        new_state.y = old_state.y + (y_af - old_state.y) * alpha_f_inv;
        new_state.ydot = old_state.ydot + (ydot_am - old_state.ydot) / alpha_m_inv;

        return new_state;
    }
} // namespace ALGEBRA

#endif // SVZERODSOLVER_ALGEBRA_INTEGRATOR_HPP_