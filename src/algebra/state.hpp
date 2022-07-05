/**
 * @file state.hpp
 * @brief ALGEBRA::State source file
 */
#ifndef SVZERODSOLVER_ALGEBRA_STATE_HPP_
#define SVZERODSOLVER_ALGEBRA_STATE_HPP_

namespace ALGEBRA
{
    /**
     * @brief State of the system.
     *
     * Stores the current state of a system, i.e. the current value and
     * derivate of all variables.
     *
     * @tparam T Scalar type (e.g. `float`, `double`)
     */
    template <typename T>
    class State
    {
    public:
        Eigen::Matrix<T, Eigen::Dynamic, 1> y;    ///< Vector of solution quantities
        Eigen::Matrix<T, Eigen::Dynamic, 1> ydot; ///< Derivate of \ref y

        /**
         * @brief Construct a new State object
         *
         */
        State();

        /**
         * @brief Construct a new State object
         *
         * @param n Size of the state
         */
        State(unsigned int n);

        /**
         * @brief Destroy the State object
         *
         */
        ~State();

        /**
         * @brief Copy a State object
         *
         * @param state
         */
        State(const State &state);

        /**
         * @brief Construct a new State object and initilaize with all zeros.
         *
         * @param n Size of the state
         * @return New state initialized with all zeros
         */
        static State Zero(unsigned int n);
    };

    template <typename T>
    State<T>::State()
    {
    }

    template <typename T>
    State<T>::State(unsigned int n)
    {
        y = Eigen::Matrix<T, Eigen::Dynamic, 1>(n);
        ydot = Eigen::Matrix<T, Eigen::Dynamic, 1>(n);
    }

    template <typename T>
    State<T>::~State()
    {
    }

    template <typename T>
    State<T>::State(const State &state)
    {
        y = state.y;
        ydot = state.ydot;
    }

    template <typename T>
    State<T> State<T>::Zero(unsigned int n)
    {
        static State<T> state(n);
        state.y = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n);
        state.ydot = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n);
        return state;
    }
} // namespace ALGEBRA

#endif // SVZERODSOLVER_ALGEBRA_STATE_HPP_