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

Integrator::Integrator(double rho, unsigned int n, double time_step_size)
{
    alpha_m = 0.5 * (3.0 - rho) / (1.0 + rho);
    alpha_f = 1.0 / (1.0 + rho);
    gamma = 0.5 + alpha_m - alpha_f;
    gamma_inv = 1.0 / gamma;

    this->n = n;
    this->time_step_size = time_step_size;
    time_step_size_inv = 1.0 / time_step_size;

    y_dot_coeff = alpha_m / (alpha_f * gamma) * time_step_size_inv;
}

Integrator::~Integrator()
{
}

State Integrator::step(State &state, double time, std::map<std::string, std::variant<Junction, BloodVessel, FlowReference, RCRBlockWithDistalPressure>> &blocks, unsigned int max_iter)
{
    State old_state = state;
    State new_state;

    new_state.y = old_state.y + 0.5 * time_step_size * old_state.ydot;
    new_state.ydot = old_state.ydot * (gamma - 0.5) * gamma_inv;

    Eigen::VectorXd y_af = old_state.y + alpha_f * (new_state.y - old_state.y);
    Eigen::VectorXd ydot_am = old_state.ydot + alpha_m * (new_state.ydot - old_state.ydot);

    double new_time = time + alpha_f * time_step_size;

    // Update time in blocks
    for (auto &&elem : blocks)
    {
        std::visit([&](auto &&block)
                   {block.update_constant(system); block.update_time(system, time); },
                   elem.second);
    }

    for (size_t i = 0; i < max_iter; i++)
    {

        // Update solution and assemble
        for (auto &&elem : blocks)
        {
            std::visit([&](auto &&block)
                       { block.update_solution(system, y_af); },
                       elem.second);
        }

        // Calculate RHS and LHS
        auto lhs = system.F + (system.dE + system.dF + system.dC + system.E * y_dot_coeff);
        auto rhs = -system.E.transpose() * ydot_am - system.F.transpose() * y_af - system.C;

        // Solve system
        Eigen::VectorXd dy = lhs.colPivHouseholderQr().solve(rhs);

        // Update solution
        y_af = y_af + dy;
        ydot_am = ydot_am * y_dot_coeff;

        if (rhs.norm() < 5.0e-4)
        {
            break;
        }
    }

    new_state.y = old_state.y + (y_af - old_state.y) / alpha_f;
    new_state.ydot = old_state.ydot + (ydot_am - old_state.ydot) / alpha_m;

    return new_state;
}