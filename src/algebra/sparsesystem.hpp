#ifndef SVZERODSOLVER_ALGREBRA_SPARSESYSTEM_HPP_
#define SVZERODSOLVER_ALGREBRA_SPARSESYSTEM_HPP_

#include <Eigen/Sparse>
#include <Eigen/SparseLU>

namespace ALGEBRA
{

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

} // namespace ALGEBRA

#endif // SVZERODSOLVER_ALGREBRA_SPARSESYSTEM_HPP_
