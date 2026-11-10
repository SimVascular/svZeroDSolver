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
#include <utility>
#include <vector>

#include "ActivationFunction.h"
#include "Block.h"
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
 * ### Constitutive law (\cite genet23 Eqs. 22-29)
 *
 * The second Piola-Kirchhoff stress is split into passive deviatoric, bulk,
 * viscous and active parts, \f$\Sigma = \Sigma^d + \Sigma^b + \Sigma^v +
 * \Sigma^a\f$, with
 * - transversely isotropic deviatoric energy
 *   \f$W^e = C_1(\bar I_1-3) + C_2(\bar I_2-3) + C_3 e^{C_4(\bar I_1-3)^2} +
 *   C_5 e^{C_6(\bar I_4-1)^2}\f$, \f$\Sigma^d=\partial W^e/\partial E\f$,
 *   optionally plus the Holzapfel-Ogden energy of svMultiPhysics
 *   (`HolzapfelOgden`)
 *   \f$W^\text{HO}=\frac{a}{2b}e^{b(\bar I_1-3)}
 *   +\sum_{i=f,s}\frac{a_i}{2b_i}\chi_i\left(e^{b_i(\bar I_{4i}-1)^2}-1\right)
 *   +\frac{a_{fs}}{2b_{fs}}\left(e^{b_{fs}\bar I_8^2}-1\right)\f$ with sheet
 *   \f$\mathbf{e}_S=\mathbf{e}_R\f$, \f$\bar I_{4f}=\bar I_4\f$,
 *   \f$\bar I_{4s}=J^{-2/3}\mathbf{e}_S\cdot C\cdot\mathbf{e}_S\f$,
 *   \f$\bar I_8=J^{-2/3}\mathbf{e}_F\cdot C\cdot\mathbf{e}_S\f$ and the
 * smoothed Heaviside \f$\chi_i=1/(1+e^{-k(\bar I_{4i}-1)})\f$;
 * - an incompressibility bulk term \f$\Sigma^b = p\,J\,C^{-1}\f$ with
 *   \f$p\f$ an independent hydrostatic-pressure field (\cite genet23 mixed
 *   u/p formulation), discretized element-wise constant (P0, LBB-stable with
 *   the linear \f$\rho,\varphi,\eta\f$) and constrained by the weak
 *   incompressibility \f$\int_{\Omega_e}(J-1)\,\mathrm d\Omega=0\f$ per element
 *   (locking-free at coarse meshes);
 * - a viscous term \f$\Sigma^v=\gamma\,\dot E\f$;
 * - an active fiber stress
 *   \f$\Sigma^a = \tau\,\mathbf{e}_F\otimes\mathbf{e}_F\f$ along the local
 *   myofiber direction \f$\mathbf{e}_F(R) = (0,\cos\alpha(R),\sin\alpha(R))\f$
 *   with helix angle varying linearly through the wall,
 *   \f$\alpha(R)=\alpha_i\frac{R_e-R}{R_e-R_i}+\alpha_e\frac{R-R_i}{R_e-R_i}\f$.
 *
 * The invariants are those of the isochoric tensor \f$\bar C = J^{-2/3}C\f$,
 * with \f$\bar I_4 = J^{-2/3}\,\mathbf{e}_F\cdot C\cdot\mathbf{e}_F\f$.
 *
 * ### Active contraction (\cite genet23 / Bestel-Clement-Sorine, Ref. 6)
 *
 * The active fiber stress is \f$\sigma_\text{1D} = (\tau_c+\mu\dot e_c)/
 * (1+e_\text{fib})\f$ (Eqs. 30-32 of the paper / Eq. 10 of Chapelle et al.
 * 2012), where \f$1+e_\text{fib}=\sqrt{\mathbf{e}_F\cdot C\cdot\mathbf{e}_F}=
 * \sqrt{I_4}\f$ from the kinematics, and the contractile-element strain
 * \f$e_c\f$, active stress \f$\tau_c\f$ and active stiffness \f$k_c\f$ evolve
 * per quadrature point via
 * \f[
 * \mu\dot e_c = k_s(e_\text{fib}-e_c)-\tau_c,\quad
 * \dot\tau_c = -(|\nu|+\alpha|\dot e_c|)\tau_c + \dot e_c\,k_c +
 *   n_0(e_c)\sigma_0|\nu|_+,\quad
 * \dot k_c = -(|\nu|+\alpha|\dot e_c|)k_c + n_0(e_c)k_0|\nu|_+,
 * \f]
 * (Chapelle Eq. 9) driven by the ECG-derived activation \f$\nu(t)\f$ (the
 * authors' piecewise-linear `nagumo` waveform, supplied as an
 * `activation_function` object; see \ref NagumoActivation),
 * \f$\sigma_0=\f$ `sigma_max`, and the piecewise-linear Frank-Starling factor
 * \f$n_0(e_c)\f$ (overridable via a `frank_starling` config object).
 * Writing the stress through \f$\tau_c+\mu\dot e_c\f$ rather than the
 * equivalent \f$k_s(e_\text{fib}-e_c)\f$ keeps the stiff \f$k_s=10^8\f$ out of
 * the mechanical residual.
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
 * \f$P_\text{in}=P_\text{out}=P_v\f$. Without inertia (`density` = 0) the rigid
 * rotation and axial translation are removed by the pins
 * \f$\varphi(R_i)=\eta(R_i)=0\f$; with inertia they are free and fixed by
 * momentum conservation (\cite genet23 Sec. 2.4).
 *
 * The ends are closed by rigid lids that carry the cavity pressure
 * (\cite genet23 Eqs. 39-44), so the wall takes the axial force
 * \f$P_v\pi(R_i+\rho(R_i))^2\f$. With `fixed_length` = 1 the ends instead slide
 * freely in fixed planes: \f$\varepsilon=0\f$ and the lids carry no load.
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
 * * `C1` ... `C6` - Passive material constants \f$C_1\ldots C_6\f$
 * * `gamma` - Material viscosity \f$\gamma\f$
 * * `sigma_max` - Maximum active fiber stress \f$\sigma_0\f$ of the BCS model
 * * `num_elements` - Number of (linear) finite elements through the wall
 *   (optional, default 10)
 * * `k_s` - BCS series-spring stiffness \f$k_s\f$ (optional, default 1e8)
 * * `k_0` - BCS maximum active stiffness \f$k_0\f$ (optional, default 260e3)
 * * `mu` - BCS active dissipation \f$\mu\f$ (optional, default 70). Distinct
 *   from the passive `gamma`.
 * * `bcs_alpha` - BCS velocity-dependent bond-destruction (force-velocity)
 *   constant \f$\alpha\f$ (optional, default 0; the paper uses 12). The
 *   \f$\alpha|\dot e_c|\f$ term is stiff and destabilizes plain
 *   generalized-alpha; set the simulation parameter `"integrator": "stiff"`
 *   (ConsistentStiffIntegrator) to use \f$\alpha>0\f$. The \f$k_c\dot e_c\f$
 *   source is retained regardless of \f$\alpha\f$.
 * * `density` - Reference mass density \f$\varrho_0\f$ (optional, default 1000
 *   kg/m^3 = the paper's 1 kg/L)
 * * `fv_reg` - BCS force-velocity regularizer \f$\epsilon\f$ smoothing
 *   \f$|\dot e_c|\approx\sqrt{\dot e_c^2+\epsilon}\f$ (optional, default 1e-4)
 * * `fixed_length` - 1 fixes the length (\f$\varepsilon=0\f$, ends sliding in
 *   fixed planes), 0 closes the ends with loaded lids (optional, default 0)
 * * `ho_a`, `ho_b`, `ho_af`, `ho_bf`, `ho_as`, `ho_bs`, `ho_afs`, `ho_bfs`,
 *   `ho_k` - Holzapfel-Ogden constants \f$a, b, a_f, b_f, a_s, b_s, a_{fs},
 *   b_{fs}, k\f$ (svMultiPhysics `a`, `b`, `a4f`, `b4f`, `a4s`, `b4s`, `afs`,
 *   `bfs`, `k`); a term is off when its \f$a\f$ is 0, otherwise its \f$b\f$
 *   must be positive (optional, default 0)
 * * `bulk_modulus` - Bulk modulus \f$\kappa\f$ of the volumetric energy
 *   \f$U=\tfrac\kappa4(J^2-1-2\ln J)\f$ (ST91, as in svMultiPhysics),
 *   enforced with the element pressure; 0 makes the wall incompressible
 *   (optional, default 0)
 *
 * The ECG-derived activation \f$\nu(t)\f$ is supplied separately as an
 * `activation_function` object (type `nagumo`, see \ref NagumoActivation) with
 * parameters `tsys`, `tdias`, `qrs`, `ramp`, `alpha_max`, `alpha_min`. An
 * optional `frank_starling` object with arrays `e_c` and `n0` overrides the
 * default force-length curve \f$n_0(e_c)\f$.
 *
 * ### Internal variables
 *
 * * `rho_<A>` - Radial displacement field nodal values \f$\rho_A\f$
 * * `phi_<A>` - In-plane shear field nodal values \f$\varphi_A\f$
 * * `eta_<A>` - Out-of-plane shear field nodal values \f$\eta_A\f$
 * * `beta` - Twist per unit length \f$\beta\f$
 * * `eps` - Longitudinal strain \f$\varepsilon\f$
 * * `ec_<q>`, `tauc_<q>`, `kc_<q>` - BCS contractile-element strain, active
 *   bond stress and active stiffness at quadrature point `q`
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
    C1,
    C2,
    C3,
    C4,
    C5,
    C6,
    gamma,
    sigma_max,
    num_elements,
    k_s,           // BCS series-spring stiffness
    k_0,           // BCS maximum active stiffness
    mu,            // BCS active dissipation
    bcs_alpha,     // BCS activation rate constant (paper's alpha)
    density,       // reference mass density rho_0 (inertia, genet23 Eq. 18)
    fv_reg,        // BCS force-velocity regularizer (|e_c dot| smoothing eps)
    fixed_length,  // 1: eps = 0 and no lid load; 0: pressure-loaded lids
    bulk_modulus,  // ST91 bulk modulus kappa; 0: incompressible
    ho_a,          // Holzapfel-Ogden isotropic term
    ho_b,
    ho_af,  // Holzapfel-Ogden fiber term
    ho_bf,
    ho_as,  // Holzapfel-Ogden sheet term
    ho_bs,
    ho_afs,  // Holzapfel-Ogden fiber-sheet term
    ho_bfs,
    ho_k  // Holzapfel-Ogden smoothed-Heaviside slope
  };

  /**
   * @brief Construct a new ChamberCylinder object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  ChamberCylinder(int id, Model* model)
      : Block(id, model, BlockType::chamber_cylinder, BlockClass::vessel,
              {{"Ri", InputParameter()},
               {"Re", InputParameter()},
               {"length", InputParameter()},
               {"alpha_endo", InputParameter()},
               {"alpha_epi", InputParameter()},
               {"C1", InputParameter()},
               {"C2", InputParameter()},
               {"C3", InputParameter()},
               {"C4", InputParameter()},
               {"C5", InputParameter()},
               {"C6", InputParameter()},
               {"gamma", InputParameter()},
               {"sigma_max", InputParameter()},
               {"num_elements", InputParameter(true, false, true, 10.0)},
               {"k_s", InputParameter(true, false, true, 1.0e8)},
               {"k_0", InputParameter(true, false, true, 260.0e3)},
               {"mu", InputParameter(true, false, true, 70.0)},
               {"bcs_alpha", InputParameter(true, false, true, 0.0)},
               {"density", InputParameter(true, false, true, 1000.0)},
               {"fv_reg", InputParameter(true, false, true, 1.0e-4)},
               {"fixed_length", InputParameter(true, false, true, 0.0)},
               {"bulk_modulus", InputParameter(true, false, true, 0.0)},
               {"ho_a", InputParameter(true, false, true, 0.0)},
               {"ho_b", InputParameter(true, false, true, 0.0)},
               {"ho_af", InputParameter(true, false, true, 0.0)},
               {"ho_bf", InputParameter(true, false, true, 0.0)},
               {"ho_as", InputParameter(true, false, true, 0.0)},
               {"ho_bs", InputParameter(true, false, true, 0.0)},
               {"ho_afs", InputParameter(true, false, true, 0.0)},
               {"ho_bfs", InputParameter(true, false, true, 0.0)},
               {"ho_k", InputParameter(true, false, true, 0.0)}}) {}

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
   * @brief Set the activation function driving the BCS active law nu(t).
   *
   * Takes ownership of the activation function (type `nagumo`), whose
   * compute(time) returns the signed activation rate nu(t).
   *
   * @param af Unique pointer to the activation function
   */
  void set_activation_function(
      std::unique_ptr<ActivationFunction> af) override {
    activation_function_ = std::move(af);
  }

  /**
   * @brief Override the Frank-Starling force-length curve n_0(e_c).
   *
   * Replaces the default piecewise-linear breakpoints. `ec` (abscissas) and
   * `val` (ordinates) must have equal length >= 2.
   *
   * @param ec Contractile-element strain breakpoints e_c
   * @param val Recruitment factor values n_0 at those breakpoints
   */
  void set_frank_starling_curve(std::vector<double> ec,
                                std::vector<double> val);

 private:
  // --- Finite element discretization of the wall thickness ---
  int n_ele = 0;   ///< Number of finite elements through the wall
  int n_node = 0;  ///< Number of nodes (n_ele + 1)
  int n_quad = 0;  ///< Number of quadrature points through the wall
  /// Rigid-body pins phi(R_i) = eta(R_i) = 0. genet23 (Sec 2.4) needs them only
  /// in the static setting; with inertia (density > 0) the rigid rotation and
  /// axial translation are free and governed by momentum conservation.
  bool pin_rigid_ = true;
  /// Fixed length (parameter `fixed_length`): eps = 0 replaces the eps
  /// equation.
  bool fixed_length_ = false;
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
  //          (e_c,tau_c,k_c)_q.., volume, vrho.., vphi.., veta.., vbeta, veps,
  //          pmix_0..]
  //   eqns: [rho.., phi.., eta.., beta, eps, (e_c,tau_c,k_c)_q.., volume, mass,
  //          pressure, vrho.., vphi.., veta.., vbeta, veps, pmix_0..]
  // The BCS block holds 3 * n_quad variables/equations (e_c, tau_c, k_c per
  // quadrature point).
  int i_rho(int a) const { return 4 + a; }
  int i_phi(int a) const { return 4 + n_node + a; }
  int i_eta(int a) const { return 4 + 2 * n_node + a; }
  int i_beta() const { return 4 + 3 * n_node; }
  int i_eps() const { return 4 + 3 * n_node + 1; }
  /// Start of the BCS variables
  int i_active0() const { return 4 + 3 * n_node + 2; }
  int i_ec(int q) const { return i_active0() + 3 * q; }
  int i_tauc(int q) const { return i_active0() + 3 * q + 1; }
  int i_kc(int q) const { return i_active0() + 3 * q + 2; }
  int n_active_var() const { return 3 * n_quad; }
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
  /// Start of the BCS equations
  int e_active0() const { return 3 * n_node + 2; }
  int e_ec(int q) const { return e_active0() + 3 * q; }
  int e_tauc(int q) const { return e_active0() + 3 * q + 1; }
  int e_kc(int q) const { return e_active0() + 3 * q + 2; }
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

  double act = 0.0;       ///< |nu(t)|
  double act_plus = 0.0;  ///< |nu(t)|_+ = max(nu(t), 0)

  /// Activation function nu(t) driving the BCS active law (type `nagumo`).
  std::unique_ptr<ActivationFunction> activation_function_;

  /// Frank-Starling breakpoints of n_0(e_c). Default: the curve of the genet23
  /// authors' input code (\cite genet23 Table 1); overridable via
  /// set_frank_starling_curve().
  std::vector<double> n0_ec_ = {-0.1668, -0.0073, 0.0534, 0.0969, 0.1326,
                                0.2016,  0.4663,  0.9187, 1.1762};
  std::vector<double> n0_val_ = {0.0, 0.5614, 0.7748, 0.8933, 0.9618,
                                 1.0, 1.0,    0.1075, 0.0};
};

#endif  // SVZERODSOLVER_MODEL_ChamberCylinder_HPP_
