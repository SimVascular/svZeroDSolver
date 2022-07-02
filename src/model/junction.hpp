/**
 * @file junction.hpp
 * @brief MODEL::Junction source file
 */
#ifndef SVZERODSOLVER_MODEL_JUNCTION_HPP_
#define SVZERODSOLVER_MODEL_JUNCTION_HPP_

#include "../algebra/densesystem.hpp"
#include "../algebra/sparsesystem.hpp"
#include "block.hpp"

namespace MODEL
{
    /**
     * @brief Junction
     *
     * Models a junction with arbitrary inlets and outlets. Across all inlets and
     * outlets of the junction, mass is conserved and pressure is continuous.
     *
     * \f[
     * \begin{circuitikz}
     * \draw node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
     * \draw (1,0) node[anchor=south]{$P_{in}$} to [short, *-*] (3.0,0);
     * \draw (3,0) node[anchor=south]{} to [short, -*] (4.5,1.0);
     * \draw (4.3,1.1) node[anchor=south] {$P_{out,1}$};
     * \draw (3,0) node[anchor=south]{} to [short, -*] (4.5,-1.0);
     * \draw (4.3,-1.1) node[anchor=north] {$P_{out,2}$};
     * \draw [-latex] (4.65,1.1) -- (5.25,1.5) node[right] {$Q_{out,1}$};
     * \draw [-latex] (4.65,-1.1) -- (5.25,-1.5) node[right] {$Q_{out,2}$};
     * \end{circuitikz}
     * \f]
     *
     * ### Governing equations
     *
     * \f[
     * \sum_{i}^{n_{inlets}} Q_{in, i}=\sum_{j}^{n_{outlets}} Q_{out, j}
     * \f]
     *
     * \f[
     * P_{i}=P_{j} \quad \mathrm{with} \quad i \neq j
     * \f]
     *
     * ### Local contributions
     *
     * \f[
     * \mathbf{y}^{e}=\left[\begin{array}{llllllllll}P_{in, 1}^{e} & Q_{in, 1}^{e} & \dots & P_{in, i}^{e} & Q_{in, i}^{e} & P_{out, 1}^{e} & Q_{out, 1}^{e} & \dots & P_{out, i}^{e} & Q_{out, i}^{e}\end{array}\right]
     * \f]
     *
     * Mass conservation
     *
     * \f[
     * \mathbf{F}^{e}_1 = \left[\begin{array}{llllllllll}0 & 1 & 0 & 1 & \dots & 0 & -1 & 0 & -1 & \dots\end{array}\right]
     *  \f]
     *
     * Due to the pressure continuity, we can write for all independent pressure pairs:
     * \f[
     * \mathbf{F}^{e}_{2,...,n} = \left[\begin{array}{lllll}\dots & \underbrace{1}_{P_i} & \dots & \underbrace{1}_{P_j} & \dots\end{array}\right] \quad \mathrm{with} \quad i \neq j
     * \f]
     *
     * @tparam T Scalar type (e.g. `float`, `double`)
     */
    template <typename T>
    class Junction : public Block<T>
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
        };

        /**
         * @brief Construct a new Junction object
         *
         * @param name Name
         */
        Junction(std::string name);

        /**
         * @brief Destroy the Junction object
         *
         */
        ~Junction();

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
         * @brief Update the constant contributions of the element in a sparse system
         *
         * @param system System to update contributions at
         */
        void update_constant(ALGEBRA::SparseSystem<T> &system);

        /**
         * @brief Number of triplets of element
         *
         * Number of triplets that the element contributes to the global system (relevant for sparse memory reservation)
         */
        std::map<std::string, int> num_triplets = {
            {"F", 0},
            {"E", 0},
            {"D", 0},
        };

    private:
        Parameters params;
        unsigned int num_inlets;
        unsigned int num_outlets;
    };

    template <typename T>
    Junction<T>::Junction(std::string name) : Block<T>(name)
    {
        this->name = name;
    }

    template <typename T>
    Junction<T>::~Junction()
    {
    }

    template <typename T>
    void Junction<T>::setup_dofs(DOFHandler &dofhandler)
    {
        // Set number of equations of a junction block based on number of
        // inlets/outlets. Must be set before calling parent constructor
        num_inlets = this->inlet_nodes.size();
        num_outlets = this->outlet_nodes.size();
        Block<T>::setup_dofs_(dofhandler, num_inlets + num_outlets, 0);
    }

    template <typename T>
    void Junction<T>::update_constant(ALGEBRA::DenseSystem<T> &system)
    {
        // Continuous pressure condition
        for (size_t i = 0; i < (num_inlets + num_outlets - 1); i++)
        {
            system.F(this->global_eqn_ids[i], this->global_var_ids[0]) = 1.0;
            system.F(this->global_eqn_ids[i], this->global_var_ids[2 * i + 2]) = -1.0;
            num_triplets["F"] += 2;
        }

        // Conservation of mass
        for (size_t i = 1; i < num_inlets * 2; i = i + 2)
        {
            system.F(this->global_eqn_ids[num_inlets + num_outlets - 1], this->global_var_ids[i]) = 1.0;
            num_triplets["F"] += 1;
        }
        for (size_t i = (num_inlets * 2) + 1; i < (num_inlets + num_outlets) * 2; i = i + 2)
        {
            system.F(this->global_eqn_ids[num_inlets + num_outlets - 1], this->global_var_ids[i]) = -1.0;
            num_triplets["F"] += 1;
        }
    }

    template <typename T>
    void Junction<T>::update_constant(ALGEBRA::SparseSystem<T> &system)
    {
        for (size_t i = 0; i < (num_inlets + num_outlets - 1); i++)
        {
            system.F.coeffRef(this->global_eqn_ids[i], this->global_var_ids[0]) = 1.0;
            system.F.coeffRef(this->global_eqn_ids[i], this->global_var_ids[2 * i + 2]) = -1.0;
        }
        for (size_t i = 1; i < num_inlets * 2; i = i + 2)
        {
            system.F.coeffRef(this->global_eqn_ids[num_inlets + num_outlets - 1], this->global_var_ids[i]) = 1.0;
        }
        for (size_t i = (num_inlets * 2) + 1; i < (num_inlets + num_outlets) * 2; i = i + 2)
        {
            system.F.coeffRef(this->global_eqn_ids[num_inlets + num_outlets - 1], this->global_var_ids[i]) = -1.0;
        }
    }

} // namespace MODEL

#endif // SVZERODSOLVER_MODEL_JUNCTION_HPP_