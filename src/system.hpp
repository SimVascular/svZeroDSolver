#ifndef SVZERODSOLVER_SYSTEM_H_
#define SVZERODSOLVER_SYSTEM_H_

#include <Eigen/Dense>

template <typename T>
class System
{
public:
    System();
    System(unsigned int n);
    ~System();
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> F;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> E;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> dF;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> dE;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> dC;
    Eigen::Matrix<T, Eigen::Dynamic, 1> C;
};

template <typename T>
System<T>::System()
{
}

template <typename T>
System<T>::System(unsigned int n)
{
    F = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(n, n);
    E = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(n, n);
    dF = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(n, n);
    dE = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(n, n);
    dC = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(n, n);
    C = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n);
}

template <typename T>
System<T>::~System()
{
}

#endif // SVZERODSOLVER_SYSTEM_H_
