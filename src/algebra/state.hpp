#ifndef SVZERODSOLVER_ALGEBRA_STATE_HPP_
#define SVZERODSOLVER_ALGEBRA_STATE_HPP_

namespace ALGEBRA
{

    template <typename T>
    class State
    {
    public:
        Eigen::Matrix<T, Eigen::Dynamic, 1> y;
        Eigen::Matrix<T, Eigen::Dynamic, 1> ydot;
        State();
        ~State();
        State(const State &state);

        static State Zero(unsigned int n);
    };
} // namespace ALGEBRA

#endif // SVZERODSOLVER_ALGEBRA_STATE_HPP_