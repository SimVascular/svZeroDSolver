#include "system.hpp"

System::System()
{
}

System::~System()
{
}

void System::setup_matrices(const unsigned int n)
{
    F = Eigen::MatrixXd::Zero(n, n);
    E = Eigen::MatrixXd::Zero(n, n);
    dF = Eigen::MatrixXd::Zero(n, n);
    dE = Eigen::MatrixXd::Zero(n, n);
    dC = Eigen::MatrixXd::Zero(n, n);
    C = Eigen::VectorXd::Zero(n);
}