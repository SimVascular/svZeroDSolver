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
 * Models the mechanical behavior of a spherical heart chamber with active contraction.
 *
 * ### Helper Functions
 * 
 * Cauchy-Green deformation tensor and time derivative:
 * \f[
 * C(r) = (1 + \frac{r}{r_0})^2
 * \dot{C}(r, \dot{r}) = 2 * (1 + \frac{r}{r_0}) \frac{\dot{r}}{r_0}
 * \f]
 *
 * ### Governing equations
 *
 * 1. Balance of linear momentum:
 * \f[
 * \rho d_0 \dot{v} + (d_0 / r_0) (1 + \frac{r}{r_0}) S - P_\text{out} C(r) = 0
 * \f]
 *
 * 2. Spherical stress:
 * \f[
 * -S + \tau + 4 (1 - C(r)^{-3}) (W_1 + C(r) W_2) + 2 \eta \dot{C}(r, \dot{r}) (1 - 2 C(r)^{-6}) = 0
 * \f]
 *
 * 3. Volume change:
 * \f[
 * 4 \pi r_0^2 C(r)v - \dot{V} = 0
 * \f]
 *
 * 4. Active stress:
 * \f[
 * \dot{\tau} + a \tau - sigma_\text{max} a_+ = 0, \quad a_+ = \max(a, 0)
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
 * \text{Pin} - \text{Pout} = 0
 * \f]
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block:
 *
 * * `rho` - Density
 * * `thick0` - Wall thickness
 * * `radius0` - Reference radius
 * * `W1` - Material constant 1
 * * `W2` - Material constant 2
 * * `eta` - Viscosity parameter
 * * `sigma_max` - Maximum active stress
 * * `alpha_max` - Maximum activation parameter
 * * `alpha_min` - Minimum activation parameter
 * * `tsys` - Systole timing parameter
 * * `tdias` - Diastole timing parameter
 * * `steepness` - Activation steepness parameter
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
