/**
 * @file openloopcoronarybc.hpp
 * @brief MODEL::ClosedLoopCoronaryBC source file
 */
#ifndef SVZERODSOLVER_MODEL_CLOSEDLOOPCORONARYBC_HPP_
#define SVZERODSOLVER_MODEL_CLOSEDLOOPCORONARYBC_HPP_

#include "../algebra/densesystem.hpp"
#include "../algebra/sparsesystem.hpp"
#include "block.hpp"

namespace MODEL {

/**
 * @brief Closed loop coronary boundary condition (connected to other blocks on both sides).
 *
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0) node[anchor=south]{$P_{in}$}
 * to [R, l=$R_a$, *-] (3,0)
 * to [R, l=$R_{am}$, -] (5,0)
 * to [R, l=$R_v$, *-*] (7,0)
 * node[anchor=south]{$P_{out}$}
 * node[left] {$Q_{out}$} [-latex] (7.2,0) -- (8.0,0);
 * (5,0) to [C, l=$C_{im} \;V_{im}$, -*] (5,-1.5)
 * node[left]{$P_{im}$}
 * (3,0) to [C, l=$C_a$, -*] (3,-1.5)
 * node[left]{$P_a$};
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f[
 * C_{i m} R_{v} Q^{e}-V_{i m}^{e}-C_{i m} P_{i m}+C_{i m} P_{v}-C_{i m} R_{v}
 * \frac{d V_{i m}^{e}}{d t}-C_{a} C_{i m} R_{v} \frac{d P^{e}}{d t}+R_{a} C_{a}
 * C_{i m} R_{v} \frac{d Q^{e}}{d t}+C_{a} C_{i m} R_{v} \frac{d P_{a}^{e}}{d
 * t}=0 \f]
 *
 * \f[
 * C_{i m} R_v P^{e}-C_{i m} R_{v} R_{a} Q^{e}-R_{v} V_{i m}^{e}-C_{i m} R_{v}
 * P_{i m}-C_{i m} R_{v} R_{a m} \frac{d V_{i m}^{e}}{d t}-R_{a m} V_{i
 * m}^{e}-C_{i m} R_{a m} P_{i m}+R_{a m} C_{i m} P_{v}=0 \f]
 *
 * ### Local contributions
 *
 * \f[
 * \mathbf{y}^{e}=\left[\begin{array}{lll}P^{e} & Q^{e} & V_{i
 * m}^{e}\end{array}\right]^{T}, \f]
 *
 * \f[
 * \mathbf{E}^{e}=\left[\begin{array}{ccc}-C_{a} C_{i m} R_{v} & R_{a} C_{a}
 * C_{i m} R_{v} & -C_{i m} R_{v} \\ 0 & 0 & -C_{i m} R_{v} R_{a
 * m}\end{array}\right] \f]
 *
 * \f[
 * \mathbf{F}^{e}=\left[\begin{array}{ccc}0 & C_{i m} R_{v} & -1 \\C_{i m} R_{v}
 * & -C_{i m} R_{v} R_{a} & -\left(R_{v}+R_{a m}\right)\end{array}\right] \f]
 *
 * \f[
 * \mathbf{c}^{e}=\left[\begin{array}{c}C_{i m}\left(-P_{i m}+P_{v}\right)+C_{a}
 * C_{i m} R_{v} \frac{d P_{a}}{d t} \\-C_{i m}\left(R_{v}+R_{a m}\right) P_{i
 * m}+R_{a m} C_{i m} P_{v}\end{array}\right] \f]
 *
 * Assume \f$P_a=0\f$.
 *
 * @tparam T Scalar type (e.g. `float`, `double`)
 */
template <typename T>
class ClosedLoopCoronaryBC : public Block<T> {
 public:
  /**
   * @brief Parameters of the element.
   *
   * Struct containing all constant and/or time-dependent parameters of the
   * element.
   */
  struct Parameters : public Block<T>::Parameters {
    T Ra;
    T Ram;
    T Rv;
    T Ca;
    T Cim;
  };

  /**
   * @brief Construct a new ClosedLoopCoronaryBC object
   *
   * @param P Time dependent pressure
   * @param name Name
   */
  ClosedLoopCoronaryBC(T Ra, T Ram, T Rv, T Ca, T Cim,
                     std::string side, std::string name);

  /**
   * @brief Destroy the ClosedLoopCoronaryBC object
   *
   */
  ~ClosedLoopCoronaryBC();

  /**
   * @brief Set up the degrees of freedom (DOF) of the block
   *
   * Set \ref global_var_ids and \ref global_eqn_ids of the element based on the
   * number of equations and the number of internal variables of the
   * element.
   *
   * @param dofhandler Degree-of-freedom handler to register variables and
   * equations at
   */
  void setup_dofs(DOFHandler &dofhandler);

  /**
   * @brief Update the constant contributions of the element in a dense system
   *
   * @param system System to update contributions at
   */
  void update_constant(ALGEBRA::DenseSystem<T> &system);

  /**
   * @brief Update the solution-dependent contributions of the element
   *
   * @param system System to update contributions at
   * @param y Current solution
   */
  void update_solution(ALGEBRA::DenseSystem<T> &system,
                       Eigen::Matrix<T, Eigen::Dynamic, 1> &y);

  /**
   * @brief Update the constant contributions of the element in a sparse system
   *
   * @param system System to update contributions at
   */
  void update_constant(ALGEBRA::SparseSystem<T> &system);

  /**
   * @brief Update the solution-dependent contributions of the element
   *
   * @param system System to update contributions at
   * @param y Current solution
   */
  void update_solution(ALGEBRA::SparseSystem<T> &system,
                       Eigen::Matrix<T, Eigen::Dynamic, 1> &y);

  /**
   * @brief Update the solution-dependent contributions of t
   *
   * @param system System to update contributions at
   * @param y Current solution
   */
  void update_model_dependent_params(MODEL::Model<T> &model);

  /**
   * @brief Number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  std::map<std::string, int> num_triplets = {
      {"F", 9},
      {"E", 5},
      {"D", 0},
  };

 private:
  Parameters params;
  std::string side;      // Left or right coronary?
  int ventricle_var_id;  // Index of either left or right ventricle
  T im;                  // Either iml or imr based on left or right
};

template <typename T>
ClosedLoopCoronaryBC<T>::ClosedLoopCoronaryBC(T Ra, T Ram, T Rv, T Ca, T Cim,
                                          std::string side, std::string name)
    : Block<T>(name) {
  this->name = name;
  this->params.Ra = Ra;
  this->params.Ram = Ram;
  this->params.Rv = Rv;
  this->params.Ca = Ca;
  this->params.Cim = Cim;
  this->side = side;
  this->closed_loop_outlet = true;
}

template <typename T>
ClosedLoopCoronaryBC<T>::~ClosedLoopCoronaryBC() {}

template <typename T>
void ClosedLoopCoronaryBC<T>::setup_dofs(DOFHandler &dofhandler) {
  Block<T>::setup_dofs_(dofhandler, 3, 1);
}

template <typename T>
void ClosedLoopCoronaryBC<T>::update_constant(ALGEBRA::DenseSystem<T> &system) {
  system.E(this->global_eqn_ids[0], this->global_var_ids[0]) = -params.Ram*params.Ca;
  system.E(this->global_eqn_ids[0], this->global_var_ids[1]) = params.Ram*params.Ra*params.Ca;
  system.E(this->global_eqn_ids[1], this->global_var_ids[0]) = -params.Ca;
  system.E(this->global_eqn_ids[1], this->global_var_ids[1]) = params.Ca*params.Ra;
  system.E(this->global_eqn_ids[1], this->global_var_ids[4]) = -1.0;
  
  system.F(this->global_eqn_ids[0], this->global_var_ids[0]) = -1.0;
  system.F(this->global_eqn_ids[0], this->global_var_ids[1]) = (params.Ra + params.Ram);
  system.F(this->global_eqn_ids[0], this->global_var_ids[2]) = 1.0;
  system.F(this->global_eqn_ids[0], this->global_var_ids[3]) = params.Rv;
  system.F(this->global_eqn_ids[1], this->global_var_ids[1]) = 1.0;
  system.F(this->global_eqn_ids[1], this->global_var_ids[3]) = -1.0;
  system.F(this->global_eqn_ids[2], this->global_var_ids[2]) = params.Cim;
  system.F(this->global_eqn_ids[2], this->global_var_ids[3]) = params.Cim*params.Rv;
  system.F(this->global_eqn_ids[2], this->global_var_ids[4]) = -1.0;
}

template <typename T>
void ClosedLoopCoronaryBC<T>::update_solution(ALGEBRA::DenseSystem<T> &system,
                                     Eigen::Matrix<T, Eigen::Dynamic, 1> &y) {
  auto Pim = this->im*y[this->ventricle_var_id];
  system.C(this->global_eqn_ids[2]) = -params.Cim*Pim;
}

template <typename T>
void ClosedLoopCoronaryBC<T>::update_constant(ALGEBRA::SparseSystem<T> &system) {
  system.E.coeffRef(this->global_eqn_ids[0], this->global_var_ids[0]) = -params.Ram*params.Ca;
  system.E.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) = params.Ram*params.Ra*params.Ca;
  system.E.coeffRef(this->global_eqn_ids[1], this->global_var_ids[0]) = -params.Ca;
  system.E.coeffRef(this->global_eqn_ids[1], this->global_var_ids[1]) = params.Ca*params.Ra;
  system.E.coeffRef(this->global_eqn_ids[1], this->global_var_ids[4]) = -1.0;
  
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[0]) = -1.0;
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) = (params.Ra + params.Ram);
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[2]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[3]) = params.Rv;
  system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[1]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[3]) = -1.0;
  system.F.coeffRef(this->global_eqn_ids[2], this->global_var_ids[2]) = params.Cim;
  system.F.coeffRef(this->global_eqn_ids[2], this->global_var_ids[3]) = params.Cim*params.Rv;
  system.F.coeffRef(this->global_eqn_ids[2], this->global_var_ids[4]) = -1.0;
}

template <typename T>
void ClosedLoopCoronaryBC<T>::update_solution(ALGEBRA::SparseSystem<T> &system,
                                     Eigen::Matrix<T, Eigen::Dynamic, 1> &y) {
  auto Pim = this->im*y[this->ventricle_var_id];
  system.C(this->global_eqn_ids[2]) = -params.Cim*Pim;
}

template <typename T>
void ClosedLoopCoronaryBC<T>::update_model_dependent_params(MODEL::Model<T> &model) {
  T im_value = 0.0;
  for (auto &[key, elem] : model.blocks) {
    std::visit([&](auto &&block) { 
      if (key == "CLH0") {
        if (this->side == "left") {
          block.get_parameter_value("iml", im_value); // Scaling for LV pressure -> intramyocardial pressure
          this->ventricle_var_id = block.global_var_ids[13]; // Solution ID for LV pressure
        }
        else if (this->side == "right") {
          block.get_parameter_value("imr", im_value); // Scaling for RV pressure -> intramyocardial pressure
          this->ventricle_var_id = block.global_var_ids[6]; 
        }
        else {
          throw std::runtime_error("For closed loop coronary, 'side' should be either 'left' or 'right'");
        }
      }
    }, elem);
  }
  this->im = im_value;
}

}  // namespace MODEL

#endif  // SVZERODSOLVER_MODEL_OPENLOOPCORONARYBC_HPP_
