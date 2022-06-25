#ifndef SVZERODSOLVER_SYSTEM_H_
#define SVZERODSOLVER_SYSTEM_H_

#include <Eigen/Dense>

class System
{
public:
    System();
    System(unsigned int n);
    ~System();
    Eigen::MatrixXd F;
    Eigen::MatrixXd E;
    Eigen::MatrixXd dF;
    Eigen::MatrixXd dE;
    Eigen::MatrixXd dC;
    Eigen::VectorXd C;
};

#endif // SVZERODSOLVER_SYSTEM_H_
