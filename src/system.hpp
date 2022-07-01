#ifndef SVZERODSOLVER_SYSTEM_H_
#define SVZERODSOLVER_SYSTEM_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

template <typename T>
class System
{
public:
    System();
    System(unsigned int n);
    ~System();
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> F;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> E;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> D;
    Eigen::Matrix<T, Eigen::Dynamic, 1> C;

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> jacobian;
    Eigen::Matrix<T, Eigen::Dynamic, 1> residual;
    Eigen::Matrix<T, Eigen::Dynamic, 1> dy;

    void reserve(std::map<std::string, int> num_triplets);
    void update_residual(Eigen::Matrix<T, Eigen::Dynamic, 1> &y, Eigen::Matrix<T, Eigen::Dynamic, 1> &ydot);
    void update_jacobian(T e_coeff);
    void solve();
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
    D = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(n, n);
    C = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n);

    jacobian = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(n, n);
    residual = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n);
    dy = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n);
}

template <typename T>
System<T>::~System()
{
}

template <typename T>
void System<T>::reserve(std::map<std::string, int> num_triplets)
{
}

template <typename T>
void System<T>::update_residual(Eigen::Matrix<T, Eigen::Dynamic, 1> &y, Eigen::Matrix<T, Eigen::Dynamic, 1> &ydot)
{
    residual = -(E * ydot) - (F * y) - C;
}

template <typename T>
void System<T>::update_jacobian(T e_coeff)
{
    jacobian = F + D + E * e_coeff;
}

template <typename T>
void System<T>::solve()
{
    // TODO: Works only if matrix is invertable: Check if True otherwise use colPivHouseholderQr
    dy = jacobian.partialPivLu().solve(residual);
}

template <typename T>
class SparseSystem
{
public:
    SparseSystem();
    SparseSystem(unsigned int n);
    ~SparseSystem();
    Eigen::SparseMatrix<T> F;
    Eigen::SparseMatrix<T> E;
    Eigen::SparseMatrix<T> D;
    Eigen::Matrix<T, Eigen::Dynamic, 1> C;

    Eigen::SparseMatrix<T> jacobian;
    Eigen::Matrix<T, Eigen::Dynamic, 1> residual;
    Eigen::Matrix<T, Eigen::Dynamic, 1> dy;

    void reserve(std::map<std::string, int> num_triplets);
    void update_residual(Eigen::Matrix<T, Eigen::Dynamic, 1> &y, Eigen::Matrix<T, Eigen::Dynamic, 1> &ydot);
    void update_jacobian(T e_coeff);
    void solve();
};

template <typename T>
SparseSystem<T>::SparseSystem()
{
}

template <typename T>
SparseSystem<T>::SparseSystem(unsigned int n)
{
    F = Eigen::SparseMatrix<T>(n, n);
    E = Eigen::SparseMatrix<T>(n, n);
    D = Eigen::SparseMatrix<T>(n, n);
    C = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n);

    jacobian = Eigen::SparseMatrix<T>(n, n);
    residual = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n);
    dy = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n);
}

template <typename T>
SparseSystem<T>::~SparseSystem()
{
}

template <typename T>
void SparseSystem<T>::reserve(std::map<std::string, int> num_triplets)
{
    F.reserve(num_triplets["F"]);
    E.reserve(num_triplets["E"]);
    D.reserve(num_triplets["D"]);
    jacobian.reserve(num_triplets["F"] + num_triplets["E"]);
}

template <typename T>
void SparseSystem<T>::update_residual(Eigen::Matrix<T, Eigen::Dynamic, 1> &y, Eigen::Matrix<T, Eigen::Dynamic, 1> &ydot)
{
    residual = -(E * ydot) - (F * y) - C;
}

template <typename T>
void SparseSystem<T>::update_jacobian(T e_coeff)
{
    jacobian = (F + D) + (E * e_coeff);
}

template <typename T>
void SparseSystem<T>::solve()
{
    Eigen::SparseLU<Eigen::SparseMatrix<T>> solver;
    solver.compute(jacobian);
    dy = solver.solve(residual);
}

#endif // SVZERODSOLVER_SYSTEM_H_
