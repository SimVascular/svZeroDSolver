// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file ChamberCylinder.h
 * @brief model::ChamberCylinder source file
 */
#ifndef SVZERODSOLVER_MODEL_ChamberCylinder_HPP_
#define SVZERODSOLVER_MODEL_ChamberCylinder_HPP_

#include <math.h>

#include <Eigen/Dense>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "ActivationFunction.h"
#include "ActiveStress.h"
#include "Block.h"
#include "Material.h"
#include "SparseSystem.h"

/**
 * @brief Thick-walled cylindrical heart chamber model
 *
 * Reduced-order left-ventricular model based on the cylindrical geometry and
 * kinematics of \cite genet23 ("Reduced left ventricular dynamics modeling
 * based on a cylindrical assumption"). The ventricle is represented as a
 * closed thick-walled cylinder of reference internal radius \f$R_i\f$, external
 * radius \f$R_e\f$ and length \f$L\f$. In the spirit of the spherical
 * ChamberSphere block, this block only models the ventricular wall mechanics
 * and exposes a single cavity pressure at its inlet/outlet; valves and the
 * circulation are modeled by separate 0D blocks.
 *
 * Compared to ChamberSphere, this model resolves the myofiber orientation and
 * the mechanical state **through the wall thickness** by integrating the
 * governing equations over the radial direction \f$R\in[R_i,R_e]\f$ with
 * **finite elements**, and it captures ventricular **twist**.
 *
 * ### Kinematics (\cite genet23 Eq. 1)
 *
 * The admissible deformation is
 * \f[
 * r = R + \rho(R),\quad \theta = \Theta + \beta Z + \varphi(R),\quad
 * z = (1+\varepsilon) Z + \eta(R),
 * \f]
 * parametrized by the three radial fields \f$\rho,\varphi,\eta\f$ (discretized
 * with Lagrange finite elements across the thickness) and the two global
 * scalars \f$\beta\f$ (twist per unit length) and \f$\varepsilon\f$
 * (longitudinal strain). The deformation gradient, right Cauchy-Green tensor
 * \f$C\f$ and Green-Lagrange strain \f$E=\tfrac12(C-1)\f$ follow \cite genet23
 * Eqs. (13,15,16), with volume ratio
 * \f$J=(1+\rho')(1+\rho/R)(1+\varepsilon)\f$.
 *
 * ### Constitutive law (\cite genet23 Eqs. 22-33)
 *
 * The second Piola-Kirchhoff stress is split into passive deviatoric, bulk,
 * viscous and active parts, \f$\Sigma = \Sigma^d + \Sigma^b + \Sigma^v +
 * \Sigma^a\f$, with
 * - a deviatoric stress \f$\Sigma^d=2\,\partial W/\partial C\f$ of the wall
 *   \ref Material (`material` object), with the fiber direction
 *   \f$\mathbf{e}_F(R) = (0,\cos\alpha(R),\sin\alpha(R))\f$ and the radial
 *   sheet direction \f$\mathbf{e}_S=\mathbf{e}_R\f$ in the
 *   \f$(\mathbf{e}_R,\mathbf{e}_\Theta,\mathbf{e}_Z)\f$ basis; the helix
 *   angle varies linearly through the wall,
 *   \f$\alpha(R)=\alpha_i\frac{R_e-R}{R_e-R_i}+\alpha_e\frac{R-R_i}{R_e-R_i}\f$;
 * - a bulk term \f$\Sigma^b = p\,J\,C^{-1}\f$ with \f$p\f$ an independent
 *   hydrostatic-pressure field (\cite genet23 mixed u/p formulation),
 *   discretized element-wise constant (P0, LBB-stable with the linear
 *   \f$\rho,\varphi,\eta\f$) and determined by the incompressibility
 *   constraint \f$\int_{\Omega_e}(J-1)\,\mathrm d\Omega=0\f$ of each
 *   element (locking-free on coarse meshes);
 * - a viscous term \f$\Sigma^v=\gamma\,\dot E\f$;
 * - an active fiber stress
 *   \f$\Sigma^a = \tau\,\mathbf{e}_F\otimes\mathbf{e}_F\f$ of the
 *   \ref ActiveStress law (`active_stress` object), evaluated at each
 *   quadrature point at the fiber stretch
 *   \f$\lambda=\sqrt{\mathbf{e}_F\cdot C\cdot\mathbf{e}_F}\f$ and driven by
 *   the `activation_function` (unless the law takes none).
 *
 * The material and active stress laws are those of \ref ChamberSphere. The
 * genet23 model is the material `mooney_rivlin_exponential` (\cite genet23
 * Eq. 25) with the active stress `genet23` (the Bestel-Clement-Sorine law of
 * Eqs. 30-33) driven by an activation rate (`piecewise_rate`).
 *
 * ### Governing equations
 *
 * The dynamic principle of virtual work (\cite genet23 Eqs. 8, 18, 45), reduced
 * to one spatial dimension (\cite genet23 Appendix A.2, Eqs. A10-A12), reads
 * \f[
 * P_a + 2\pi L\int_{R_i}^{R_e}\Sigma:\mathrm{D}E(\hat\zeta)\,\mathrm{d}R
 * = P_v\,2\pi L (R_i+\rho(R_i))(1+\varepsilon)\,\hat\rho(R_i)
 * + P_v\,\pi L (R_i+\rho(R_i))^2\,\hat\varepsilon \quad\forall\hat\zeta,
 * \f]
 * assembled with finite elements over the thickness, with the acceleration
 * virtual power
 * \f$P_a=\int_\Omega\varrho_0\,\ddot u\cdot\mathrm{D}u(\hat\zeta)\,d\Omega\f$
 * realized by velocity companion DOFs \f$w=\dot\zeta\f$ and the consistent
 * mass matrix
 * \f$M=\int\varrho_0(\mathrm{D}u)^\top\mathrm{D}u\,\mathrm d\Omega\f$
 * (analytically integrated over \f$\Theta,Z\f$; \cite genet23 Eqs. A5-A6).
 * It is complemented by the cavity volume
 * \f$V=\pi(R_i+\rho(R_i))^2(1+\varepsilon)L\f$, mass conservation
 * \f$Q_\text{in}-Q_\text{out}-\dot V = 0\f$ and pressure equality
 * \f$P_\text{in}=P_\text{out}=P_v\f$. The \f$\hat\varepsilon\f$ term is the
 * axial load of the rigid lids that close the tube (\cite genet23 Eqs.
 * 39-44): the wall carries \f$P_v\pi(R_i+\rho(R_i))^2\f$ and the length is
 * free. Without inertia (`density` = 0) the rigid rotation and axial
 * translation are removed by the pins \f$\varphi(R_i)=\eta(R_i)=0\f$; with
 * inertia they are free and fixed by momentum conservation (\cite genet23
 * Sec. 2.4).
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block:
 *
 * * `Ri` - Reference internal radius \f$R_i\f$
 * * `Re` - Reference external radius \f$R_e\f$
 * * `length` - Reference length \f$L\f$
 * * `alpha_endo` - Myofiber helix angle at endocardium \f$\alpha_i\f$ [deg]
 * * `alpha_epi` - Myofiber helix angle at epicardium \f$\alpha_e\f$ [deg]
 * * `gamma` - Material viscosity \f$\gamma\f$
 * * `num_elements` - Number of (linear) finite elements through the wall
 *   (optional, default 10)
 * * `density` - Reference mass density \f$\varrho_0\f$ (optional, default 1000
 *   kg/m^3 = the paper's 1 kg/L; 0 for the quasi-static model)
 *
 * Like \ref ChamberSphere, the block is declared in `chambers` with these
 * parameters in `values`, a `material` object (see \ref Material), an
 * `active_stress` object (see \ref ActiveStress) and an
 * `activation_function` object (see \ref ActivationFunction). With the
 * active stress `genet23` and \f$\alpha>0\f$ (its velocity-dependent bond
 * destruction, 12 in the paper), the stiff \f$\alpha|\dot e_c|\f$ term needs
 * the simulation parameter `"integrator": "stiff"`
 * (ConsistentStiffIntegrator).
 *
 * ### Internal variables
 *
 * * `rho_<A>` - Radial displacement field nodal values \f$\rho_A\f$
 * * `phi_<A>` - In-plane shear field nodal values \f$\varphi_A\f$
 * * `eta_<A>` - Out-of-plane shear field nodal values \f$\eta_A\f$
 * * `beta` - Twist per unit length \f$\beta\f$
 * * `eps` - Longitudinal strain \f$\varepsilon\f$
 * * `active_stress_<q>`, ... - States of the active stress law at quadrature
 *   point `q` (e.g. `active_stress_<q>`, `e_c_<q>`, `tau_c_<q>`, `k_c_<q>` for
 *   `genet23`)
 * * `volume` - Cavity volume \f$V\f$
 * * `vrho_<A>`, `vphi_<A>`, `veta_<A>`, `vbeta`, `veps` - velocity companion
 *   DOFs \f$w=\dot\zeta\f$ (full dynamics)
 * * `pmix_<e>` - Element hydrostatic pressure \f$p_e\f$ (P0 mixed u/p)
 *
 * (node index `A` runs from 0 at the endocardium \f$R_i\f$ to `num_elements`
 * at the epicardium \f$R_e\f$.)
 */
class ChamberCylinder : public Block {
 public:
  /**
   * @brief Local IDs of the parameters
   */
  enum ParamId {
    Ri,
    Re,
    length,
    alpha_endo,
    alpha_epi,
    gamma,
    num_elements,
    density  // reference mass density rho_0 (inertia, genet23 Eq. 18)
  };

  /**
   * @brief Construct a new ChamberCylinder object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  ChamberCylinder(int id, Model* model)
      : Block(id, model, BlockType::chamber_cylinder, BlockClass::chamber,
              {{"Ri", InputParameter()},
               {"Re", InputParameter()},
               {"length", InputParameter()},
               {"alpha_endo", InputParameter()},
               {"alpha_epi", InputParameter()},
               {"gamma", InputParameter()},
               {"num_elements", InputParameter(true, false, true, 10.0)},
               {"density", InputParameter(true, false, true, 1000.0)}}) {}

  /**
   * @brief Set up the degrees of freedom (DOF) of the block
   *
   * @param dofhandler Degree-of-freedom handler to register variables and
   * equations at
   */
  void setup_dofs(DOFHandler& dofhandler) override;

  /**
   * @brief Update the constant contributions of the element in a sparse system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  void update_constant(SparseSystem& system,
                       std::vector<double>& parameters) override;

  /**
   * @brief Update the time-dependent contributions of the element in a sparse
   * system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  void update_time(SparseSystem& system,
                   std::vector<double>& parameters) override;

  /**
   * @brief Update the solution-dependent contributions of the element in a
   * sparse system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   * @param y Current solution
   * @param dy Current derivate of the solution
   */
  void update_solution(
      SparseSystem& system, std::vector<double>& parameters,
      const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
      const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) override;

  /**
   * @brief Number of triplets of element
   */
  TripletsContributions num_triplets;

  /**
   * @brief Get number of triplets of element
   *
   * @return TripletsContributions Number of triplets of element
   */
  TripletsContributions get_num_triplets() override { return num_triplets; }

  /**
   * @brief Set the activation function driving the active stress law.
   *
   * @param af Unique pointer to the activation function
   */
  void set_activation_function(
      std::unique_ptr<ActivationFunction> af) override {
    activation_function_ = std::move(af);
  }

  /**
   * @brief Set the wall material.
   *
   * @param m Unique pointer to the material
   */
  void set_material(std::unique_ptr<Material> m) override {
    material_ = std::move(m);
  }

  /**
   * @brief Set the active stress law.
   *
   * @param as Unique pointer to the active stress law
   */
  void set_active_stress(std::unique_ptr<ActiveStress> as) override {
    active_stress_ = std::move(as);
  }

 private:
  // --- Finite element discretization of the wall thickness ---
  int n_ele = 0;   ///< Number of finite elements through the wall
  int n_node = 0;  ///< Number of nodes (n_ele + 1)
  int n_quad = 0;  ///< Number of quadrature points through the wall
  /// Rigid-body pins phi(R_i) = eta(R_i) = 0. genet23 (Sec 2.4) needs them only
  /// in the static setting; with inertia (density > 0) the rigid rotation and
  /// axial translation are free and governed by momentum conservation.
  bool pin_rigid_ = true;
  /// Check that the material, active stress law and activation are set.
  void check_inputs() const;
  int n_law = 0;  ///< Number of active stress states per quadrature point
  int n_var = 0;  ///< Number of block variables (external + internal)
  int n_eqn = 0;  ///< Number of block equations

  /// Reference nodal radii (endocardium at index 0, epicardium at n_node-1)
  std::vector<double> node_R;

  /**
   * @brief A single Gauss quadrature point in the thickness
   */
  struct QuadPoint {
    int elem;      ///< Element index (spans nodes elem and elem+1)
    double R;      ///< Reference radius of the quadrature point
    double w;      ///< Integration weight (Gauss weight * dR/dxi * R)
    double alpha;  ///< Myofiber helix angle at this radius [rad]
    double N[2];   ///< Shape function values at the two element nodes
    double dN[2];  ///< Shape function radial derivatives
  };
  std::vector<QuadPoint> quad;  ///< Quadrature points across the thickness

  // --- Local variable / equation index helpers ---
  // Layout:
  //   vars: [Pin, Qin, Pout, Qout, rho_0.., phi_0.., eta_0.., beta, eps,
  //          (active stress states)_q.., volume, vrho.., vphi.., veta.., vbeta,
  //          veps, pmix_0..]
  //   eqns: [rho.., phi.., eta.., beta, eps, (active stress equations)_q..,
  //          volume, mass, pressure, vrho.., vphi.., veta.., vbeta, veps,
  //          pmix_0..]
  // The active stress block holds n_law * n_quad variables/equations.
  int i_rho(int a) const { return 4 + a; }
  int i_phi(int a) const { return 4 + n_node + a; }
  int i_eta(int a) const { return 4 + 2 * n_node + a; }
  int i_beta() const { return 4 + 3 * n_node; }
  int i_eps() const { return 4 + 3 * n_node + 1; }
  /// Start of the active stress variables
  int i_active0() const { return 4 + 3 * n_node + 2; }
  /// State i of the active stress law at quadrature point q (state 0: tau)
  int i_law(int q, int i) const { return i_active0() + n_law * q + i; }
  int n_active_var() const { return n_law * n_quad; }
  int i_vol() const { return i_active0() + n_active_var(); }
  // Velocity companion DOFs (dynamics): w = d(field)/dt, appended after the
  // volume DOF.
  int i_vel0() const { return i_vol() + 1; }
  int i_vrho(int a) const { return i_vel0() + a; }
  int i_vphi(int a) const { return i_vel0() + n_node + a; }
  int i_veta(int a) const { return i_vel0() + 2 * n_node + a; }
  int i_vbeta() const { return i_vel0() + 3 * n_node; }
  int i_veps() const { return i_vel0() + 3 * n_node + 1; }
  // Mixed u/p element pressure DOFs (P0), appended after the velocity block.
  int i_pmix0() const { return i_vel0() + 3 * n_node + 2; }
  int i_pmix(int e) const { return i_pmix0() + e; }

  int e_rho(int a) const { return a; }
  int e_phi(int a) const { return n_node + a; }
  int e_eta(int a) const { return 2 * n_node + a; }
  int e_beta() const { return 3 * n_node; }
  int e_eps() const { return 3 * n_node + 1; }
  /// Start of the active stress equations
  int e_active0() const { return 3 * n_node + 2; }
  /// Equation i of the active stress law at quadrature point q
  int e_law(int q, int i) const { return e_active0() + n_law * q + i; }
  int e_volume() const { return e_active0() + n_active_var(); }
  int e_mass() const { return e_volume() + 1; }
  int e_pressure() const { return e_volume() + 2; }
  // Velocity companion equations (dynamics): w - d(field)/dt = 0.
  int e_vel0() const { return e_pressure() + 1; }
  int e_vrho(int a) const { return e_vel0() + a; }
  int e_vphi(int a) const { return e_vel0() + n_node + a; }
  int e_veta(int a) const { return e_vel0() + 2 * n_node + a; }
  int e_vbeta() const { return e_vel0() + 3 * n_node; }
  int e_veps() const { return e_vel0() + 3 * n_node + 1; }
  // Mixed u/p element incompressibility equations, appended last.
  int e_pmix0() const { return e_vel0() + 3 * n_node + 2; }
  int e_pmix(int e) const { return e_pmix0() + e; }

  // --- Local tangents and their sparsity pattern ---
  // update_solution accumulates dC/dy and dC/dydot in dense local arrays
  // (n_eqn x n_var, zero between calls) and writes only the entries its
  // assembly touches to the system. The first assembly (SparseSystem::reserve)
  // records them; the pattern is fixed afterwards, since the solver analyzes
  // it once. Which entries are touched depends on the coupling of the
  // unknowns only, not on their values.
  std::vector<double> K_;          ///< local dC/dy
  std::vector<double> Kd_;         ///< local dC/dydot
  std::vector<char> in_K_;         ///< entry of K_ in the pattern
  std::vector<char> in_Kd_;        ///< entry of Kd_ in the pattern
  std::vector<int> pattern_K_;     ///< pattern of K_ (local linear indices)
  std::vector<int> pattern_Kd_;    ///< pattern of Kd_
  bool pattern_recorded_ = false;  ///< pattern complete

  /// Wall material
  std::unique_ptr<Material> material_;

  /// Active stress law
  std::unique_ptr<ActiveStress> active_stress_;

  /// Activation function driving the active stress law
  std::unique_ptr<ActivationFunction> activation_function_;
};

#endif  // SVZERODSOLVER_MODEL_ChamberCylinder_HPP_
