/**
 * @file sparsesystem.hpp
 * @brief ALGEBRA::SparseSystem source file
 */
#ifndef SVZERODSOLVER_ALGREBRA_SPARSESYSTEM_HPP_
#define SVZERODSOLVER_ALGREBRA_SPARSESYSTEM_HPP_

#include <Eigen/Sparse>
#include <Eigen/SparseLU>

namespace ALGEBRA
{
    /**
     * @brief Sparse system
     *
     * This class contains all attributes and methods to create, modify, and
     * solve sparse systems.
     *
     * @tparam T Scalar type (e.g. `float`, `double`)
     */
    template <typename T>
    class SparseSystem
    {
    public:
        /**
         * @brief Construct a new Sparse System object
         *
         */
        SparseSystem();

        /**
         * @brief Construct a new Sparse System object
         *
         * @param n Size of the system
         */
        SparseSystem(unsigned int n);

        /**
         * @brief Destroy the Sparse System object
         *
         */
        ~SparseSystem();

        Eigen::SparseMatrix<T> F;              ///< System matrix F
        Eigen::SparseMatrix<T> E;              ///< System matrix E
        Eigen::SparseMatrix<T> D;              ///< System matrix D
        Eigen::Matrix<T, Eigen::Dynamic, 1> C; ///< System vector C

        Eigen::SparseMatrix<T> jacobian;              ///< Jacobian of the system
        Eigen::Matrix<T, Eigen::Dynamic, 1> residual; ///< Residual of the system
        Eigen::Matrix<T, Eigen::Dynamic, 1> dy;       ///< Solution increment of the system

        /**
         * @brief Reserve memory in system matrices based on number of triplets
         *
         * @param num_triplets Number of triplets that will be assembled to each matrix.
         */
        void reserve(std::map<std::string, int> num_triplets);

        /**
         * @brief Update the residual of the system
         *
         * @param y Vector of current solution quantities
         * @param ydot Derivate of y
         */
        void update_residual(Eigen::Matrix<T, Eigen::Dynamic, 1> &y, Eigen::Matrix<T, Eigen::Dynamic, 1> &ydot);

        /**
         * @brief Update the jacobian of the system
         *
         * @param e_coeff Coefficent for system matrix \ref E
         */
        void update_jacobian(T e_coeff);

        /**
         * @brief Solve the system
         */
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
