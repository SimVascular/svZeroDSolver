// Copyright (c) Stanford University, The Regents of the University of
//               California, and others.
//
// All Rights Reserved.
//
// See Copyright-SimVascular.txt for additional details.
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject
// to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
// OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/**
 * @file ChamberSphere.h
 * @brief model::ChamberSphere source file
 */
#ifndef SVZERODSOLVER_MODEL_ChamberSphere_HPP_
#define SVZERODSOLVER_MODEL_ChamberSphere_HPP_

#include <math.h>

#include "Block.h"
#include "SparseSystem.h"

/**
 * @brief Spherical heart chamber model
 *
 * Models the mechanical behavior of a spherical heart chamber with active
 * contraction. For reference, see \cite caruel13 Equations (13a-b) for
 * continuum mechanics (without length-dependent contraction valves, vessels)
 and \cite
 * pfaller2019importance Equations (12-16) for the simplified active contraction
 * model.
 *
 * ### Helper Functions
 *
 * Cauchy-Green deformation tensor and time derivative:
 * \f[
 * C = \left(1 + \frac{r}{r_0} \right)^2
 * \f]
 * \f[
 * \dot{C} = 2 \left(1 + \frac{r}{r_0} \right) \frac{\dot{r}}{r_0}
 * \f]
 *
 * ### Governing equations
 *
 * 1. Balance of linear momentum:
 * \f[
 * \rho d_0 \dot{v} + \frac{d_0}{r_0} \left(1 + \frac{r}{r_0} \right) S -
 P_\text{out} C = 0
 * \f]
 *
 * 2. Spherical stress:
 * \f[
 * -S + \tau + 4 (1 - C^{-3}) (W_1 + C W_2) + 2 \eta \dot{C}
 * (1 - 2 C^{-6}) = 0
 * \f]
 *
 * 3. Volume change:
 * \f[
 * 4 \pi r_0^2 Cv - \dot{V} = 0
 * \f]
 *
 * 4. Active stress:
 * \f[
 * \dot{\tau} + a \tau - \sigma_\text{max} a_+ = 0, \quad a_+ = \max(a, 0),
 \quad a = f\alpha_\text{max} + (1 - f)\alpha_\text{min}
 * \f]
 * with indicator function
 * \f[
 * f = S_+ \cdot S_-, \quad S_\pm = \frac{1}{2} \left(1.0 \pm \text{tanh}\left(
 \frac{t - t_\text{sys/dias}} {\gamma} \right) \right)
 * \f]
 *
 * 5. Acceleration:
 * \f[
 * \dot{r} - v = 0
 * \f]
 *
 * 6. Conservation of mass:
 * \f[
 * Q_\text{in} - Q_\text{out} - \dot{V} = 0
 * \f]
 *
 * 7. Pressure equality:
 * \f[
 * P_\text{in} - P_\text{out} = 0
 * \f]
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block:
 *
 * * `rho` - Density \f$\rho\f$
 * * `thick0` - Wall thickness \f$d_0\f$
 * * `radius0` - Reference radius \f$r_0\f$
 * * `W1` - Material constant \f$W_1\f$
 * * `W2` - Material constant \f$W_2\f$
 * * `eta` - Viscosity parameter \f$\eta\f$
 * * `sigma_max` - Maximum active stress \f$\sigma_\text{max}\f$
 * * `alpha_max` - Maximum activation parameter \f$\alpha_\text{max}\f$
 * * `alpha_min` - Minimum activation parameter \f$\alpha_\text{min}\f$
 * * `tsys` - Systole timing parameter \f$t_\text{sys}\f$
 * * `tdias` - Diastole timing parameter \f$t_\text{dias}\f$
 * * `steepness` - Activation steepness parameter \f$\gamma\f$
 *
 * ### Internal variables
 *
 * Names of internal variables in this block's output:
 *
 * * `radius` - Chamber radius \f$r\f$
 * * `velo` - Chamber velocity \f$\dot{r}\f$
 * * `stress` - Spherical stress \f$S\f$
 * * `tau` - Active stress \f$\tau\f$
 * * `volume` - Chamber volume \f$V\f$
 *
 */
class ChamberSphere : public Block {
 public:
  /**
   * @brief Local IDs of the parameters
   *
   */
  enum ParamId {
    rho = 0,
    thick0 = 1,
    radius0 = 2,
    W1 = 3,
    W2 = 4,
    eta = 5,
    sigma_max = 6,
    alpha_max = 7,
    alpha_min = 8,
    tsys = 9,
    tdias = 10,
    steepness = 11
  };

  /**
   * @brief Construct a new ChamberSphere object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  ChamberSphere(int id, Model *model)
      : Block(id, model, BlockType::chamber_sphere, BlockClass::vessel,
              {{"rho", InputParameter()},
               {"thick0", InputParameter()},
               {"radius0", InputParameter()},
               {"W1", InputParameter()},
               {"W2", InputParameter()},
               {"eta", InputParameter()},
               {"sigma_max", InputParameter()},
               {"alpha_max", InputParameter()},
               {"alpha_min", InputParameter()},
               {"tsys", InputParameter()},
               {"tdias", InputParameter()},
               {"steepness", InputParameter()}}) {}

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
   * @brief Update the constant contributions of the element in a sparse
   system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  void update_constant(SparseSystem &system, std::vector<double> &parameters);

  /**
   * @brief Update the time-dependent contributions of the element in a sparse
   * system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  void update_time(SparseSystem &system, std::vector<double> &parameters);

  /**
   * @brief Update the constant contributions of the element in a sparse
   system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  // void update_constant(SparseSystem &system, std::vector<double>
  // &parameters);

  // /**
  //  * @brief Update the solution-dependent contributions of the element in a
  //  * sparse system
  //  *
  //  * @param system System to update contributions at
  //  * @param parameters Parameters of the model
  //  * @param y Current solution
  //  * @param dy Current derivate of the solution
  //  */
  void update_solution(SparseSystem &system, std::vector<double> &parameters,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy);

  /**
   * @brief Update the elastance functions which depend on time
   *
   * @param parameters Parameters of the model
   */
  void get_elastance_values(std::vector<double> &parameters);

 private:
  double act = 0.0;       // activation function
  double act_plus = 0.0;  // act_plus = max(act, 0)

  /**
   * @brief Number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  TripletsContributions num_triplets{0, 0, 18};
};

#endif  // SVZERODSOLVER_MODEL_ChamberSphere_HPP_
