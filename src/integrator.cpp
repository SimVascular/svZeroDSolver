#include "integrator.hpp"

State::State()
{
}

State::~State()
{
}

State::State(const State &state)
{
    y = state.y;
    ydot = state.ydot;
}

State State::Zero(unsigned int n)
{
    static State state;
    state.y = Eigen::VectorXd::Zero(n);
    state.ydot = Eigen::VectorXd::Zero(n);
    return state;
}

Integrator::Integrator(System &system, double time_step_size, double rho)
{
    alpha_m = 0.5 * (3.0 - rho) / (1.0 + rho);
    alpha_f = 1.0 / (1.0 + rho);
    gamma = 0.5 + alpha_m - alpha_f;
    gamma_inv = 1.0 / gamma;

    this->system = system;
    this->time_step_size = time_step_size;
    time_step_size_inv = 1.0 / time_step_size;

    size = system.F.row(0).size();
    y_af = Eigen::VectorXd(size);
    ydot_am = Eigen::VectorXd(size);
    dy = Eigen::VectorXd(size);
    lhs = Eigen::MatrixXd(size, size);
    rhs = Eigen::VectorXd(size);

    y_dot_coeff = alpha_m / (alpha_f * gamma) * time_step_size_inv;
}

Integrator::~Integrator()
{
}

State Integrator::step(State &state, double time, Model &model, unsigned int max_iter)
{
    State old_state = state;
    State new_state;

    new_state.y = old_state.y + 0.5 * time_step_size * old_state.ydot;
    new_state.ydot = old_state.ydot * (gamma - 0.5) * gamma_inv;

    y_af = old_state.y + alpha_f * (new_state.y - old_state.y);
    ydot_am = old_state.ydot + alpha_m * (new_state.ydot - old_state.ydot);

    double new_time = time + alpha_f * time_step_size;

    // Update time in blocks
    model.update_time(system, new_time);

    dy.setZero();
    for (size_t i = 0; i < max_iter; i++)
    {
        // Update solution and assemble
        model.update_solution(system, y_af);

        // Calculate RHS and LHS
        rhs = -(system.E * ydot_am) - (system.F * y_af) - system.C;
        if (rhs.norm() < 1.0e-8)
        {
            break;
        }
        else if (i == max_iter - 1)
        {
            throw std::runtime_error("Maxium number of non-linear iterations reached.");
            exit(1);
        }
        lhs = system.F + (system.dE + system.dF + system.dC + system.E * y_dot_coeff);

        // Solve system
        dy = lhs.partialPivLu().solve(rhs);

        // Update solution
        y_af += dy;
        ydot_am += dy * y_dot_coeff;
    }

    new_state.y = old_state.y + (y_af - old_state.y) / alpha_f;
    new_state.ydot = old_state.ydot + (ydot_am - old_state.ydot) / alpha_m;

    return new_state;
}