/**
 * @file resistancebc.hpp
 * @brief MODEL::ResistanceBC source file
 */
#ifndef SVZERODSOLVER_MODEL_RESISTANCEWITHDISTALPRESSURE_HPP_
#define SVZERODSOLVER_MODEL_RESISTANCEWITHDISTALPRESSURE_HPP_

#include "../algebra/densesystem.hpp"
#include "../algebra/sparsesystem.hpp"
#include "block.hpp"
#include "timedependentparameter.hpp"

namespace MODEL
{

    /**
     * @brief Resistance boundary condition.
     *
     * \f[
     * \begin{circuitikz} \draw
     * node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
     * \draw (1.0,0) to [R, l=$R$, *-*] (3,0)
     * node[anchor=south]{$P_{d}$};
     * \end{circuitikz}
     * \f]
     *
     * ### Governing equations
     *
     * \f[
     * P-P_d=R \cdot Q
     * \f]
     *
     * ### Local contributions
     *
     * \f[
     * \mathbf{y}^{e}=\left[\begin{array}{ll}P^{e} & Q^{e}\end{array}\right]^{T}
     * \f]
     *
     * \f[
     * \mathbf{F}^{e}=\left[\begin{array}{ll}1 & -R\end{array}\right]
     * \f]
     *
     * \f[
     * \mathbf{C}^{e}=\left[-P_d\right]
     * \f]
     *
     *
     * @tparam T Scalar type (e.g. `float`, `double`)
     */
    template <typename T>
    class ResistanceBC : public Block<T>
    {
    public:
        /**
         * @brief Parameters of the element.
         *
         * Struct containing all constant and/or time-dependent parameters of the
         * element.
         */
        struct Parameters : public Block<T>::Parameters
        {
            TimeDependentParameter<T> R;  ///< Time-dependent resistance
            TimeDependentParameter<T> Pd; ///< Time-dependent distal pressure
        };

        /**
         * @brief Construct a new ResistanceBC object
         *
         * @param R Time-dependent resistance
         * @param Pd Time-dependent distal pressure
         * @param name Name
         */
        ResistanceBC(TimeDependentParameter<T> R, TimeDependentParameter<T> Pd, std::string name);

        /**
         * @brief Destroy the ResistanceBC object
         *
         */
        ~ResistanceBC();

        /**
         * @brief Set up the degrees of freedom (DOF) of the block
         *
         * Set \ref global_var_ids and \ref global_eqn_ids of the element based on the
         * number of equations and the number of internal variables of the
         * element.
         *
         * @param dofhandler Degree-of-freedom handler to register variables and equations at
         */
        void setup_dofs(DOFHandler &dofhandler);

        /**
         * @brief Update the constant contributions of the element in a dense system
         *
         * @param system System to update contributions at
         */
        void update_constant(ALGEBRA::DenseSystem<T> &system);

        /**
         * @brief Update the time-dependent contributions of the element in a dense system
         *
         * @param system System to update contributions at
         * @param time Current time
         */
        void update_time(ALGEBRA::DenseSystem<T> &system, T time);

        /**
         * @brief Update the constant contributions of the element in a sparse system
         *
         * @param system System to update contributions at
         */
        void update_constant(ALGEBRA::SparseSystem<T> &system);

        /**
         * @brief Update the time-dependent contributions of the element in a sparse system
         *
         * @param system System to update contributions at
         * @param time Current time
         */
        void update_time(ALGEBRA::SparseSystem<T> &system, T time);

        /**
         * @brief Number of triplets of element
         *
         * Number of triplets that the element contributes to the global system (relevant for sparse memory reservation)
         */
        std::map<std::string, int> num_triplets = {
            {"F", 1},
            {"E", 0},
            {"D", 0},
        };

        /**
         * @brief Convert the block to a steady behavior
         *
         * Converts the resistance and distal pressure to the constant means of
         * themselve
         */
        void to_steady();

    private:
        Parameters params;
    };

    template <typename T>
    ResistanceBC<T>::ResistanceBC(TimeDependentParameter<T> R, TimeDependentParameter<T> Pd, std::string name) : Block<T>(name)
    {
        this->name = name;
        this->params.R = R;
        this->params.Pd = Pd;
    }

    template <typename T>
    ResistanceBC<T>::~ResistanceBC()
    {
    }

    template <typename T>
    void ResistanceBC<T>::setup_dofs(DOFHandler &dofhandler)
    {
        Block<T>::setup_dofs_(dofhandler, 1, 0);
    }

    template <typename T>
    void ResistanceBC<T>::update_constant(ALGEBRA::DenseSystem<T> &system)
    {
        system.F(this->global_eqn_ids[0], this->global_var_ids[0]) = 1.0;
    }

    template <typename T>
    void ResistanceBC<T>::update_time(ALGEBRA::DenseSystem<T> &system, T time)
    {
        system.F(this->global_eqn_ids[0], this->global_var_ids[1]) = -params.R.get(time);
        system.C(this->global_eqn_ids[0]) = -params.Pd.get(time);
    }

    template <typename T>
    void ResistanceBC<T>::update_constant(ALGEBRA::SparseSystem<T> &system)
    {
        system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[0]) = 1.0;
    }

    template <typename T>
    void ResistanceBC<T>::update_time(ALGEBRA::SparseSystem<T> &system, T time)
    {
        system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) = -params.R.get(time);
        system.C(this->global_eqn_ids[0]) = -params.Pd.get(time);
    }

    template <typename T>
    void ResistanceBC<T>::to_steady()
    {
        params.R.to_steady();
        params.Pd.to_steady();
    }

} // namespace MODEL

#endif // SVZERODSOLVER_MODEL_RESISTANCEWITHDISTALPRESSURE_HPP_