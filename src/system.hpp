#ifndef SVZERODSOLVER_SYSTEM_H_
#define SVZERODSOLVER_SYSTEM_H_

#include <Eigen/Dense>

class System
{
public:
    System();
    ~System();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> F;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> E;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dF;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dE;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dC;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> C;

    void setup_matrices(unsigned int n);
};

#endif // SVZERODSOLVER_SYSTEM_H_
