/**
 * @file windkesselbc.hpp
 * @brief MODEL::WindkesselBC source file
 */
#ifndef SVZERODSOLVER_MODEL_WINDKESSELBC_HPP_
#define SVZERODSOLVER_MODEL_WINDKESSELBC_HPP_

#include "../algebra/densesystem.hpp"
#include "../algebra/sparsesystem.hpp"
#include "block.hpp"

namespace MODEL
{
    /**
     * @brief Windkessel RCR boundary condition.
     *
     * Models the mechanical behavior of a Windkessel boundary condition.
     *
     * \f[
     * \begin{circuitikz} \draw
     * node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
     * \draw (1,0) node[anchor=south]{$P_{in}$}
     * to [R, l=$R_p$, *-] (3,0)
     * to [R, l=$R_d$, *-*] (5,0)
     * node[anchor=south]{$P_{ref}$}
     * (3,0) to [C, l=$C$, *-] (3,-1.5)
     * node[ground]{$P_{C}$};
     * \end{circuitikz}
     * \f]
     *
     * ### Governing equations
     *
     * \f[
     * R_{d} Q^{e}-P_{c}^{e}+P_{r e f}-R_{d} C \frac{d P_{c}^{e}}{d t}=0
     * \f]
     *
     * \f[
     * P^{e}-P_{c}^{e}-R_{p} Q^{e}=0
     * \f]
     *
     * ### Local contributions
     *
     * \f[
     * \mathbf{y}^{e}=\left[\begin{array}{lll}P^{e} & Q^{e} & P_{c}^{e}\end{array}\right]^{T}
     * \f]
     *
     * \f[
     * \mathbf{E}^{e}=\left[\begin{array}{ccc}
     * 0 & 0 & -R_{d} C \\
     * 0 & 0 & 0
     * \end{array}\right]
     * \f]
     *
     * \f[
     * \mathbf{F}^{e}=\left[\begin{array}{ccc}
     * 0 & R_{d} & -1 \\
     * 1 & -R_{p} & -1
     * \end{array}\right]
     * \f]
     *
     * \f[
     * \mathbf{c}^{e}=\left[\begin{array}{c}
     * P_{r e f} \\
     * 0
     * \end{array}\right]
     * \f]
     *
     * @tparam T Scalar type (e.g. `float`, `double`)
     */
    template <typename T>
    class WindkesselBC : public Block<T>
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
            T Rp; ///< Proximal resistance
            T C;  ///< Capacitance
            T Rd; ///< Distal restistance
            T Pd; ///< Distal Pressure
        };

        /**
         * @brief Construct a new WindkesselBC object
         *
         * @param Rp Proximal resistance
         * @param C Capacitance
         * @param Rd Distal resistance
         * @param Pd Distal pressure
         * @param name Name
         */
        WindkesselBC(T Rp, T C, T Rd, T Pd, std::string name);

        /**
         * @brief Destroy the WindkesselBC object
         *
         */
        ~WindkesselBC();

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
            {"F", 5},
            {"E", 1},
            {"D", 0},
        };

        /**
         * @brief Convert the block to a steady behavior
         *
         * Set the capacitance to 0.
         */
        void to_steady();

    private:
        Parameters params;
    };

    template <typename T>
    WindkesselBC<T>::WindkesselBC(T Rp, T C, T Rd, T Pd, std::string name) : Block<T>(name)
    {
        this->name = name;
        this->params.Rp = Rp;
        this->params.C = C;
        this->params.Rd = Rd;
        this->params.Pd = Pd;
    }

    template <typename T>
    WindkesselBC<T>::~WindkesselBC()
    {
    }

    template <typename T>
    void WindkesselBC<T>::setup_dofs(DOFHandler &dofhandler)
    {
        Block<T>::setup_dofs_(dofhandler, 2, 1);
    }

    template <typename T>
    void WindkesselBC<T>::update_constant(ALGEBRA::DenseSystem<T> &system)
    {
        system.F(this->global_eqn_ids[0], this->global_var_ids[0]) = 1.0;
        system.F(this->global_eqn_ids[0], this->global_var_ids[2]) = -1.0;
        system.F(this->global_eqn_ids[1], this->global_var_ids[2]) = -1.0;
    }

    template <typename T>
    void WindkesselBC<T>::update_time(ALGEBRA::DenseSystem<T> &system, T time)
    {
        system.E(this->global_eqn_ids[1], this->global_var_ids[2]) = -params.Rd * params.C;
        system.F(this->global_eqn_ids[0], this->global_var_ids[1]) = -params.Rp;
        system.F(this->global_eqn_ids[1], this->global_var_ids[1]) = params.Rd;
        system.C(this->global_eqn_ids[1]) = params.Pd;
    }

    template <typename T>
    void WindkesselBC<T>::update_constant(ALGEBRA::SparseSystem<T> &system)
    {
        system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[0]) = 1.0;
        system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[2]) = -1.0;
        system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[2]) = -1.0;
    }

    template <typename T>
    void WindkesselBC<T>::update_time(ALGEBRA::SparseSystem<T> &system, T time)
    {
        system.E.coeffRef(this->global_eqn_ids[1], this->global_var_ids[2]) = -params.Rd * params.C;
        system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) = -params.Rp;
        system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[1]) = params.Rd;
        system.C(this->global_eqn_ids[1]) = params.Pd;
    }

    template <typename T>
    void WindkesselBC<T>::to_steady()
    {
        params.C = 0.0;
    }

} // namespace MODEL

#endif // SVZERODSOLVER_MODEL_WINDKESSELBC_HPP_