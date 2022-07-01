#ifndef SVZERODSOLVER_ALGREBRA_DENSESYSTEM_HPP_
#define SVZERODSOLVER_ALGREBRA_DENSESYSTEM_HPP_

#include "../external/eigen/Eigen/Dense"

namespace ALGEBRA
{

    template <typename T>
    class DenseSystem
    {
    public:
        DenseSystem();
        DenseSystem(unsigned int n);
        ~DenseSystem();
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
    DenseSystem<T>::DenseSystem()
    {
    }

    template <typename T>
    DenseSystem<T>::DenseSystem(unsigned int n)
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
    DenseSystem<T>::~DenseSystem()
    {
    }

    template <typename T>
    void DenseSystem<T>::reserve(std::map<std::string, int> num_triplets)
    {
    }

    template <typename T>
    void DenseSystem<T>::update_residual(Eigen::Matrix<T, Eigen::Dynamic, 1> &y, Eigen::Matrix<T, Eigen::Dynamic, 1> &ydot)
    {
        residual = -(E * ydot) - (F * y) - C;
    }

    template <typename T>
    void DenseSystem<T>::update_jacobian(T e_coeff)
    {
        jacobian = F + D + E * e_coeff;
    }

    template <typename T>
    void DenseSystem<T>::solve()
    {
        // TODO: Works only if matrix is invertable: Check if True otherwise use colPivHouseholderQr
        dy = jacobian.partialPivLu().solve(residual);
    }

} // namespace ALGEBRA

#endif // SVZERODSOLVER_ALGREBRA_DENSESYSTEM_HPP_
