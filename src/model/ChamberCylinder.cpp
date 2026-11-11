// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ChamberCylinder.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Model.h"

namespace {

using Eigen::Matrix3d;
using Eigen::Vector3d;

// The continuum kernel is templated on the scalar type: double for the
// residual, std::complex for complex-step tangents.
template <typename T>
using Mat3 = Eigen::Matrix<T, 3, 3>;
template <typename T>
using Vec3 = Eigen::Matrix<T, 3, 1>;
using cplx = std::complex<double>;

/// Real part of a scalar, for comparisons on complex-step arguments.
inline double re(double x) { return x; }
inline double re(const cplx& x) { return x.real(); }

/// Upper clamp compatible with complex-step: x below m, the constant m above.
template <typename T>
inline T clamp_max(T x, double m) {
  return re(x) < m ? x : T(m);
}

/// Holzapfel-Ogden constants (svMultiPhysics `HolzapfelOgden`).
struct HolzapfelOgden {
  double a = 0.0, b = 0.0;      ///< isotropic term
  double af = 0.0, bf = 0.0;    ///< fiber term
  double as = 0.0, bs = 0.0;    ///< sheet term
  double afs = 0.0, bfs = 0.0;  ///< fiber-sheet term
  double k = 0.0;               ///< smoothed-Heaviside slope
  bool active() const {
    return a != 0.0 || af != 0.0 || as != 0.0 || afs != 0.0;
  }
};

/// Material parameters passed to the continuum kernel.
struct MatParams {
  double C1, C2, C3, C4, C5, C6, gamma;
  HolzapfelOgden ho;
  double p_mix = 0.0;  ///< mixed u/p element hydrostatic pressure
};

/// Active-stress input at an integration point, giving the 2nd Piola-Kirchhoff
/// fiber stress ((tauc + mu * ecdot) / sqrt(I4) + tension) e_F (x) e_F.
struct ActiveInput {
  double tauc = 0.0;     ///< BCS active bond stress tau_c
  double ecdot = 0.0;    ///< BCS contractile-element strain rate
  double mu = 0.0;       ///< BCS active dissipation
  double tension = 0.0;  ///< prescribed active tension T(t)
};

/// Double contraction A:B (without complex conjugation, for complex-step).
template <typename T>
inline T contract(const Mat3<T>& A, const Mat3<T>& B) {
  return (A.array() * B.array()).sum();
}

/**
 * @brief Kinematics at a through-wall integration point.
 *
 * Holds the right Cauchy-Green tensor C and the six derivatives of the
 * Green-Lagrange strain with respect to the local kinematic quantities
 * xi = [rho, rho', phi', beta, eps, eta'] (\cite genet23 Eqs. 15, 16, 21).
 */
template <typename T>
struct Kinematics {
  Mat3<T> C;      ///< Right Cauchy-Green deformation tensor
  Mat3<T> dE[6];  ///< dE/dxi_k for xi=[rho,rho',phi',beta,eps,eta']
  T J;            ///< Volume ratio det(F)
};

/// Build C, J and the strain derivatives from local kinematics at radius R.
template <typename T>
Kinematics<T> compute_kinematics(const T xi[6], double R) {
  const T rho = xi[0], drho = xi[1], dphi = xi[2];
  const T beta = xi[3], eps = xi[4], deta = xi[5];
  const T a = 1.0 + drho;     // radial stretch F_RR
  const T b = 1.0 + rho / R;  // circumferential stretch F_ThTh
  const T e = 1.0 + eps;      // longitudinal stretch F_ZZ

  Kinematics<T> k;
  // Right Cauchy-Green tensor in the orthonormal (e_R, e_Th, e_Z) basis.
  const T C00 = a * a + (R * b * dphi) * (R * b * dphi) + deta * deta;
  const T C11 = b * b;
  const T C22 = (R * b * beta) * (R * b * beta) + e * e;
  const T C01 = R * b * b * dphi;
  const T C02 = R * R * b * b * dphi * beta + deta * e;
  const T C12 = R * b * b * beta;
  k.C << C00, C01, C02, C01, C11, C12, C02, C12, C22;
  k.J = a * b * e;

  for (int i = 0; i < 6; i++) k.dE[i].setZero();

  // dE/drho (with d b / d rho = 1/R)
  {
    Mat3<T>& m = k.dE[0];
    m(0, 0) = R * b * dphi * dphi;
    m(1, 1) = b / R;
    m(2, 2) = R * b * beta * beta;
    m(0, 1) = m(1, 0) = b * dphi;
    m(0, 2) = m(2, 0) = R * b * dphi * beta;
    m(1, 2) = m(2, 1) = b * beta;
  }
  // dE/drho'
  k.dE[1](0, 0) = a;
  // dE/dphi'
  {
    Mat3<T>& m = k.dE[2];
    m(0, 0) = R * R * b * b * dphi;
    m(0, 1) = m(1, 0) = 0.5 * R * b * b;
    m(0, 2) = m(2, 0) = 0.5 * R * R * b * b * beta;
  }
  // dE/dbeta
  {
    Mat3<T>& m = k.dE[3];
    m(2, 2) = R * R * b * b * beta;
    m(0, 2) = m(2, 0) = 0.5 * R * R * b * b * dphi;
    m(1, 2) = m(2, 1) = 0.5 * R * b * b;
  }
  // dE/deps
  {
    Mat3<T>& m = k.dE[4];
    m(2, 2) = e;
    m(0, 2) = m(2, 0) = 0.5 * deta;
  }
  // dE/deta'
  {
    Mat3<T>& m = k.dE[5];
    m(0, 0) = deta;
    m(0, 2) = m(2, 0) = 0.5 * e;
  }
  return k;
}

/// Fiber structural tensor e_F (x) e_F with e_F = (0, cos alpha, sin alpha) in
/// the (e_R, e_Th, e_Z) basis (\cite genet23 Eq. 22).
template <typename T>
inline Mat3<T> fiber_tensor(double alpha) {
  const double ca = std::cos(alpha), sa = std::sin(alpha);
  Vec3<T> ef(T(0.0), T(ca), T(sa));
  return ef * ef.transpose();
}

/// Symmetric fiber-sheet tensor (e_F (x) e_S + e_S (x) e_F) / 2 with the
/// sheet e_S = e_R.
template <typename T>
inline Mat3<T> fiber_sheet_tensor(double alpha) {
  Mat3<T> H = Mat3<T>::Zero();
  H(0, 1) = H(1, 0) = T(0.5 * std::cos(alpha));
  H(0, 2) = H(2, 0) = T(0.5 * std::sin(alpha));
  return H;
}

/// dW/dI of a Holzapfel-Ogden fiber or sheet energy
/// W = a/(2b) chi (exp(b E^2) - 1) with E = I - 1 and the smoothed Heaviside
/// chi = 1/(1 + exp(-k E)).
template <typename T>
inline T ho_family(T E, double a, double b, double k) {
  if (a == 0.0) return T(0.0);
  const T o = 1.0 / (std::exp(clamp_max(k * E, 300.0)) + 1.0);  // 1 - chi
  const T rexp = std::exp(clamp_max(b * E * E, 300.0));
  return a * ((1.0 - o) * E * rexp + 0.5 * k * (o - o * o) / b * (rexp - 1.0));
}

/**
 * @brief Second Piola-Kirchhoff stress at an integration point.
 *
 * Sigma = Sigma^d (deviatoric, \cite genet23 Eq. 25) + Sigma^b (hydrostatic
 * pressure, Eqs. 26-27) + Sigma^v (viscous, Eq. 28) + Sigma^a (active fiber
 * stress, Eqs. 29-30).
 */
template <typename T>
Mat3<T> compute_stress(const Kinematics<T>& k, const T xidot[6],
                       const ActiveInput& act, double alpha,
                       const MatParams& p) {
  const Mat3<T>& C = k.C;
  const T J = k.J;
  const Mat3<T> Cinv = C.inverse();
  const Mat3<T> I = Mat3<T>::Identity();
  const Mat3<T> M = fiber_tensor<T>(alpha);

  const T I1 = C.trace();
  const T trC2 = (C.array() * C.array()).sum();
  const T I2 = 0.5 * (I1 * I1 - trC2);
  const T I4 = (C * M).trace();  // e_F . C . e_F

  const T Jm23 = std::pow(J, -2.0 / 3.0);
  const T Jm43 = Jm23 * Jm23;
  const T I1b = Jm23 * I1;
  const T I2b = Jm43 * I2;
  const T I4b = Jm23 * I4;
  (void)I2b;  // I2b enters only through its derivative

  // dW/d(isochoric invariant); the exponents are capped so that non-physical
  // trial states stay finite.
  const T a1 = p.C4 * (I1b - 3.0) * (I1b - 3.0);
  const T a4 = p.C6 * (I4b - 1.0) * (I4b - 1.0);
  const T w1 =
      p.C1 + 2.0 * p.C3 * p.C4 * (I1b - 3.0) * std::exp(clamp_max(a1, 300.0));
  const T w2 = T(p.C2);
  const T w4 = 2.0 * p.C5 * p.C6 * (I4b - 1.0) * std::exp(clamp_max(a4, 300.0));

  // d(isochoric invariant)/dC
  const Mat3<T> dI1b = Jm23 * (I - (I1 / 3.0) * Cinv);
  const Mat3<T> dI2b = Jm43 * (I1 * I - C - (2.0 / 3.0) * I2 * Cinv);
  const Mat3<T> dI4b = Jm23 * (M - (I4 / 3.0) * Cinv);

  Mat3<T> Sd = 2.0 * (w1 * dI1b + w2 * dI2b + w4 * dI4b);  // deviatoric

  // Holzapfel-Ogden terms with the sheet e_S = e_R, using the isochoric
  // invariants I6b = J^-2/3 C_RR and I8b = J^-2/3 e_F . C . e_S.
  if (p.ho.active()) {
    const HolzapfelOgden& h = p.ho;
    Mat3<T> Hs = Mat3<T>::Zero();
    Hs(0, 0) = T(1.0);
    const Mat3<T> Hfs = fiber_sheet_tensor<T>(alpha);
    const T I6 = C(0, 0);
    const T I8 = contract<T>(C, Hfs);
    const T I6b = Jm23 * I6;
    const T I8b = Jm23 * I8;
    const Mat3<T> dI6b = Jm23 * (Hs - (I6 / 3.0) * Cinv);
    const Mat3<T> dI8b = Jm23 * (Hfs - (I8 / 3.0) * Cinv);
    const T wa = 0.5 * h.a * std::exp(clamp_max(h.b * (I1b - 3.0), 300.0));
    const T wf = ho_family<T>(I4b - 1.0, h.af, h.bf, h.k);
    const T ws = ho_family<T>(I6b - 1.0, h.as, h.bs, h.k);
    const T wfs = h.afs * I8b * std::exp(clamp_max(h.bfs * I8b * I8b, 300.0));
    Sd += 2.0 * (wa * dI1b + wf * dI4b + ws * dI6b + wfs * dI8b);
  }
  // Hydrostatic stress p J C^-1 with the element pressure p (mixed u/p).
  const T bulk = T(p.p_mix);
  Mat3<T> Sb = bulk * J * Cinv;  // bulk

  Mat3<T> Edot = Mat3<T>::Zero();  // viscous
  for (int m = 0; m < 6; m++) Edot += k.dE[m] * xidot[m];
  Mat3<T> Sv = p.gamma * Edot;

  // Active fiber stress sigma_1D = (tau_c + mu e_c_dot) / (1 + e_fib) with
  // 1 + e_fib = sqrt(I4) (Eqs. 30, 32). Using tau_c + mu e_c_dot instead of
  // k_s (e_fib - e_c) keeps the stiff k_s out of the mechanical residual.
  // A prescribed tension T enters the 2nd Piola-Kirchhoff stress directly,
  // S_a = T e_F (x) e_F (svMultiPhysics active stress).
  const T sigma_act =
      (act.tauc + act.mu * act.ecdot) / std::sqrt(I4) + act.tension;
  Mat3<T> Sa = sigma_act * M;

  return Sd + Sb + Sv + Sa;
}

/// Generalized internal forces g_k = Sigma : dE/dxi_k at an integration point.
template <typename T>
void compute_g(const T xi[6], const T xidot[6], const ActiveInput& act,
               double R, double alpha, const MatParams& p, T g[6]) {
  Kinematics<T> k = compute_kinematics<T>(xi, R);
  Mat3<T> S = compute_stress<T>(k, xidot, act, alpha, p);
  for (int i = 0; i < 6; i++) g[i] = contract<T>(S, k.dE[i]);
}

/// Frank-Starling factor n_0(e_c) (Caruel et al. 2013, Ref. 13 of
/// \cite genet23): linear interpolation of the breakpoints X, Y, constant
/// outside them. e_c may carry a complex-step perturbation.
template <typename T>
inline T frank_starling(T e_c, const std::vector<double>& X,
                        const std::vector<double>& Y) {
  const int NP = static_cast<int>(X.size());
  const double x = re(e_c);
  if (x <= X[0]) return T(Y[0]);
  if (x >= X[NP - 1]) return T(Y[NP - 1]);
  int k = 0;
  while (k < NP - 1 && x >= X[k + 1]) ++k;
  const double slope = (Y[k + 1] - Y[k]) / (X[k + 1] - X[k]);
  // Keep e_c in the linear term so that complex-step recovers the slope.
  return Y[k] + slope * (e_c - X[k]);
}

/// BCS active-law parameters.
struct BcsParams {
  double k_s, k_0, mu, alpha, sigma0;
  /// Regularization eps in |e_c_dot| ~ sqrt(e_c_dot^2 + eps)
  double fv_reg = 1e-4;
};

/// Residuals of the BCS ODEs (\cite genet23 Eqs. 32-33) at a quadrature point
/// for the state a = [e_c, tau_c, k_c] and its rate ad. n0_x, n0_y are the
/// Frank-Starling breakpoints.
template <typename T>
void bcs_residual(const T a[], const T ad[], T e_fib, double nu_abs,
                  double nu_plus, const BcsParams& bp,
                  const std::vector<double>& n0_x,
                  const std::vector<double>& n0_y, T res[]) {
  const T ecdot = ad[0];
  // |e_c_dot|, smoothed to be differentiable at zero
  const T abs_ecd = std::sqrt(ecdot * ecdot + bp.fv_reg);
  const T n0 = frank_starling<T>(a[0], n0_x, n0_y);
  const T decay = nu_abs + bp.alpha * abs_ecd;
  // series-spring force balance: mu * e_c_dot = k_s (e_fib - e_c) - tau_c
  res[0] = ecdot - (bp.k_s * (e_fib - a[0]) - a[1]) / bp.mu;
  // tau_c dot = -decay*tau_c + e_c_dot*k_c + n0*sigma0*nu_+
  res[1] = ad[1] - (-decay * a[1] + ecdot * a[2] + n0 * bp.sigma0 * nu_plus);
  // k_c dot = -decay*k_c + n0*k_0*nu_+
  res[2] = ad[2] - (-decay * a[2] + n0 * bp.k_0 * nu_plus);
}

}  // namespace

void ChamberCylinder::setup_dofs(DOFHandler& dofhandler) {
  // Geometry and discretization
  const double ri = model->get_parameter_value(global_param_ids[ParamId::Ri]);
  const double re = model->get_parameter_value(global_param_ids[ParamId::Re]);
  const double a_endo =
      model->get_parameter_value(global_param_ids[ParamId::alpha_endo]);
  const double a_epi =
      model->get_parameter_value(global_param_ids[ParamId::alpha_epi]);
  n_ele = std::max(1, (int)std::lround(model->get_parameter_value(
                          global_param_ids[ParamId::num_elements])));
  n_node = n_ele + 1;
  pin_rigid_ =
      model->get_parameter_value(global_param_ids[ParamId::density]) == 0.0;
  fixed_length_ =
      model->get_parameter_value(global_param_ids[ParamId::fixed_length]) > 0.5;

  // Uniform nodes from endocardium (index 0, R_i) to epicardium (R_e).
  node_R.resize(n_node);
  const double h = (re - ri) / n_ele;
  for (int i = 0; i < n_node; i++) node_R[i] = ri + i * h;

  // Three-point Gauss quadrature per element; the weight includes the radius
  // R of the volume element R dR dTheta dZ.
  const double gp[3] = {-std::sqrt(3.0 / 5.0), 0.0, std::sqrt(3.0 / 5.0)};
  const double gw[3] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
  const double deg2rad = M_PI / 180.0;
  quad.clear();
  for (int e = 0; e < n_ele; e++) {
    const double Ra = node_R[e], Rb = node_R[e + 1];
    for (int q = 0; q < 3; q++) {
      QuadPoint pt;
      pt.elem = e;
      pt.R = 0.5 * (Ra + Rb) + 0.5 * (Rb - Ra) * gp[q];
      pt.w = gw[q] * 0.5 * (Rb - Ra) * pt.R;  // Gauss weight * dR/dxi * R
      const double frac = (pt.R - ri) / (re - ri);
      pt.alpha = deg2rad * (a_endo * (1.0 - frac) + a_epi * frac);
      pt.N[0] = (Rb - pt.R) / (Rb - Ra);
      pt.N[1] = (pt.R - Ra) / (Rb - Ra);
      pt.dN[0] = -1.0 / (Rb - Ra);
      pt.dN[1] = 1.0 / (Rb - Ra);
      quad.push_back(pt);
    }
  }
  n_quad = (int)quad.size();

  // Internal variables: rho/phi/eta fields, beta, eps, active block, volume.
  std::list<std::string> int_names;
  for (int a = 0; a < n_node; a++)
    int_names.push_back("rho_" + std::to_string(a));
  for (int a = 0; a < n_node; a++)
    int_names.push_back("phi_" + std::to_string(a));
  for (int a = 0; a < n_node; a++)
    int_names.push_back("eta_" + std::to_string(a));
  int_names.push_back("beta");
  int_names.push_back("eps");
  // BCS states e_c, tau_c, k_c per quadrature point
  for (int q = 0; q < n_quad; q++) {
    int_names.push_back("ec_" + std::to_string(q));
    int_names.push_back("tauc_" + std::to_string(q));
    int_names.push_back("kc_" + std::to_string(q));
  }
  int_names.push_back("volume");

  // Velocity companions w = d(field)/dt for the inertia terms
  for (int a = 0; a < n_node; a++)
    int_names.push_back("vrho_" + std::to_string(a));
  for (int a = 0; a < n_node; a++)
    int_names.push_back("vphi_" + std::to_string(a));
  for (int a = 0; a < n_node; a++)
    int_names.push_back("veta_" + std::to_string(a));
  int_names.push_back("vbeta");
  int_names.push_back("veps");

  // Mixed u/p: one constant (P0) hydrostatic pressure per element
  for (int e = 0; e < n_ele; e++)
    int_names.push_back("pmix_" + std::to_string(e));

  // 3*n_node field eqns + beta + eps + active eqns + volume + mass + pressure
  // + velocity companion eqns + element incompressibility constraints.
  const int n_vel = 3 * n_node + 2;
  const int n_pmix = n_ele;
  n_eqn = 3 * n_node + 2 + n_active_var() + 3 + n_vel + n_pmix;
  Block::setup_dofs_(dofhandler, n_eqn, int_names);
  n_var = 4 + 3 * n_node + 2 + n_active_var() + 1 + n_vel + n_pmix;

  // Reservation estimates: dense dC_dy/dC_dydot block, few F/E entries
  num_triplets.D = n_eqn * n_var;
  num_triplets.F = 6 * n_node + 20;
  num_triplets.E = 21 * n_node + 8;  // consistent mass + velocity companions
}

void ChamberCylinder::update_constant(SparseSystem& system,
                                      std::vector<double>& parameters) {
  // Mass conservation Qin - Qout - Vdot = 0 (Vdot in update_solution)
  system.F.coeffRef(global_eqn_ids[e_mass()], global_var_ids[1]) = 1.0;
  system.F.coeffRef(global_eqn_ids[e_mass()], global_var_ids[3]) = -1.0;

  // Pressure equality Pin - Pout = 0
  system.F.coeffRef(global_eqn_ids[e_pressure()], global_var_ids[0]) = 1.0;
  system.F.coeffRef(global_eqn_ids[e_pressure()], global_var_ids[2]) = -1.0;

  // Cavity volume V - pi (R_i + rho_i)^2 (1 + eps) L = 0 (linear part)
  system.F.coeffRef(global_eqn_ids[e_volume()], global_var_ids[i_vol()]) = 1.0;

  // Rigid-body pins phi(R_i) = eta(R_i) = 0, only without inertia (genet23
  // Sec 2.4); with inertia, momentum conservation fixes the rigid modes.
  if (pin_rigid_) {
    system.F.coeffRef(global_eqn_ids[e_phi(0)], global_var_ids[i_phi(0)]) = 1.0;
    system.F.coeffRef(global_eqn_ids[e_eta(0)], global_var_ids[i_eta(0)]) = 1.0;
  }

  // Fixed length: eps = 0 replaces the eps equation.
  if (fixed_length_) {
    system.F.coeffRef(global_eqn_ids[e_eps()], global_var_ids[i_eps()]) = 1.0;
  }

  // --- Inertia (genet23 Eqs. 8, 18, A5-A6) ---
  // Companion DOFs w = zeta_dot make the system first order. The consistent
  // mass M = int rho0 Du^T Du dOmega, integrated over Theta and Z, couples
  // rho-rho, phi-phi, beta-beta, phi-beta, eta-eta, eps-eps and eta-eps. This
  // is its reference-radius part; update_solution adds the current-radius
  // correction and the convective term.
  {
    const double rho0 = parameters[global_param_ids[ParamId::density]];
    const double Lp = parameters[global_param_ids[ParamId::length]];
    // 2 pi rho0 times the Z-moments int_Z {1, Z, Z^2} = {L, L^2/2, L^3/3}
    const double cL1 = rho0 * 2.0 * M_PI * Lp;
    const double cL2 = rho0 * 2.0 * M_PI * Lp * Lp / 2.0;
    const double cL3 = rho0 * 2.0 * M_PI * Lp * Lp * Lp / 3.0;
    // Accumulate, then assign: update_constant can run more than once.
    std::map<std::pair<int, int>, double> M;
    auto Madd = [&](int r, int c, double v) {
      if (!(fixed_length_ && r == e_eps())) M[{r, c}] += v;
    };
    for (int q = 0; q < n_quad; q++) {
      const QuadPoint& pt = quad[q];
      const int A[2] = {pt.elem, pt.elem + 1};
      const double N[2] = {pt.N[0], pt.N[1]};
      const double R2 = pt.R * pt.R;
      Madd(e_beta(), i_vbeta(), cL3 * pt.w * R2);  // beta-beta
      Madd(e_eps(), i_veps(), cL3 * pt.w);         // eps-eps
      for (int i = 0; i < 2; i++) {
        // When pinned, the pins replace the phi_0/eta_0 momentum rows.
        const bool free_row = !pin_rigid_ || A[i] != 0;
        Madd(e_beta(), i_vphi(A[i]), cL2 * pt.w * R2 * N[i]);  // beta-phi
        Madd(e_eps(), i_veta(A[i]), cL2 * pt.w * N[i]);        // eps-eta
        if (free_row) {
          Madd(e_phi(A[i]), i_vbeta(), cL2 * pt.w * R2 * N[i]);
          Madd(e_eta(A[i]), i_veps(), cL2 * pt.w * N[i]);
        }
        for (int j = 0; j < 2; j++) {
          const double mrr = cL1 * pt.w * N[i] * N[j];  // rho-rho / eta-eta
          Madd(e_rho(A[i]), i_vrho(A[j]), mrr);
          if (free_row) {
            Madd(e_eta(A[i]), i_veta(A[j]), mrr);
            Madd(e_phi(A[i]), i_vphi(A[j]), cL1 * pt.w * R2 * N[i] * N[j]);
          }
        }
      }
    }
    for (const auto& kv : M)
      system.E.coeffRef(global_eqn_ids[kv.first.first],
                        global_var_ids[kv.first.second]) = kv.second;
    // Velocity companion equations w_k - d(field_k)/dt = 0
    auto companion = [&](int erow, int ivel, int ifield) {
      system.F.coeffRef(global_eqn_ids[erow], global_var_ids[ivel]) = 1.0;
      system.E.coeffRef(global_eqn_ids[erow], global_var_ids[ifield]) = -1.0;
    };
    for (int a = 0; a < n_node; a++) {
      companion(e_vrho(a), i_vrho(a), i_rho(a));
      companion(e_vphi(a), i_vphi(a), i_phi(a));
      companion(e_veta(a), i_veta(a), i_eta(a));
    }
    companion(e_vbeta(), i_vbeta(), i_beta());
    companion(e_veps(), i_veps(), i_eps());
  }
}

void ChamberCylinder::update_time(SparseSystem& system,
                                  std::vector<double>& parameters) {
  (void)parameters;
  if (active_tension_.defined()) {
    tension = active_tension_.value(model->time);
    if (activation_function_ == nullptr) {  // BCS at rest
      act = act_plus = 0.0;
      return;
    }
  }
  // BCS inputs |nu| and |nu|_+ from the signed activation rate nu(t)
  if (activation_function_ == nullptr) {
    throw std::runtime_error(
        "ChamberCylinder block '" + get_name() +
        "': no activation_function set. Provide an 'activation_function' "
        "(type 'nagumo') on the vessel configuration.");
  }
  const double nu = activation_function_->compute(model->time);
  act = std::abs(nu);
  act_plus = std::max(nu, 0.0);
}

void ChamberCylinder::update_solution(
    SparseSystem& system, std::vector<double>& parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) {
  MatParams p;
  p.C1 = parameters[global_param_ids[ParamId::C1]];
  p.C2 = parameters[global_param_ids[ParamId::C2]];
  p.C3 = parameters[global_param_ids[ParamId::C3]];
  p.C4 = parameters[global_param_ids[ParamId::C4]];
  p.C5 = parameters[global_param_ids[ParamId::C5]];
  p.C6 = parameters[global_param_ids[ParamId::C6]];
  p.gamma = parameters[global_param_ids[ParamId::gamma]];
  p.ho.a = parameters[global_param_ids[ParamId::ho_a]];
  p.ho.b = parameters[global_param_ids[ParamId::ho_b]];
  p.ho.af = parameters[global_param_ids[ParamId::ho_af]];
  p.ho.bf = parameters[global_param_ids[ParamId::ho_bf]];
  p.ho.as = parameters[global_param_ids[ParamId::ho_as]];
  p.ho.bs = parameters[global_param_ids[ParamId::ho_bs]];
  p.ho.afs = parameters[global_param_ids[ParamId::ho_afs]];
  p.ho.bfs = parameters[global_param_ids[ParamId::ho_bfs]];
  p.ho.k = parameters[global_param_ids[ParamId::ho_k]];
  const double sigma_max = parameters[global_param_ids[ParamId::sigma_max]];
  const double Rip = parameters[global_param_ids[ParamId::Ri]];
  const double Lp = parameters[global_param_ids[ParamId::length]];
  const double kappa = parameters[global_param_ids[ParamId::bulk_modulus]];

  BcsParams bp;
  bp.k_s = parameters[global_param_ids[ParamId::k_s]];
  bp.k_0 = parameters[global_param_ids[ParamId::k_0]];
  bp.mu = parameters[global_param_ids[ParamId::mu]];
  bp.alpha = parameters[global_param_ids[ParamId::bcs_alpha]];
  bp.sigma0 = sigma_max;
  bp.fv_reg = parameters[global_param_ids[ParamId::fv_reg]];

  // Block DOFs and their rates
  auto Y = [&](int lv) { return y[global_var_ids[lv]]; };
  auto Yd = [&](int lv) { return dy[global_var_ids[lv]]; };
  const double beta = Y(i_beta()), eps = Y(i_eps());
  const double beta_d = Yd(i_beta()), eps_d = Yd(i_eps());
  const double Pin = Y(0);

  // Local dense residual and tangents, scattered into the system at the end
  std::vector<double> Cloc(n_eqn, 0.0);
  std::vector<double> Kloc(n_eqn * n_var, 0.0);
  std::vector<double> Kdloc(n_eqn * n_var, 0.0);
  auto Kat = [&](int r, int c) -> double& { return Kloc[r * n_var + c]; };
  auto Kdat = [&](int r, int c) -> double& { return Kdloc[r * n_var + c]; };

  // --- Finite element assembly of the internal virtual work over thickness ---
  for (int q = 0; q < n_quad; q++) {
    const QuadPoint& pt = quad[q];
    const int A0 = pt.elem, A1 = pt.elem + 1;
    const double N0 = pt.N[0], N1 = pt.N[1], dN0 = pt.dN[0], dN1 = pt.dN[1];

    // Local kinematics xi = [rho, rho', phi', beta, eps, eta'] and rates
    const double rho = N0 * Y(i_rho(A0)) + N1 * Y(i_rho(A1));
    const double drho = dN0 * Y(i_rho(A0)) + dN1 * Y(i_rho(A1));
    const double dphi = dN0 * Y(i_phi(A0)) + dN1 * Y(i_phi(A1));
    const double deta = dN0 * Y(i_eta(A0)) + dN1 * Y(i_eta(A1));
    const double xi[6] = {rho, drho, dphi, beta, eps, deta};

    const double rho_d = N0 * Yd(i_rho(A0)) + N1 * Yd(i_rho(A1));
    const double drho_d = dN0 * Yd(i_rho(A0)) + dN1 * Yd(i_rho(A1));
    const double dphi_d = dN0 * Yd(i_phi(A0)) + dN1 * Yd(i_phi(A1));
    const double deta_d = dN0 * Yd(i_eta(A0)) + dN1 * Yd(i_eta(A1));
    const double xidot[6] = {rho_d, drho_d, dphi_d, beta_d, eps_d, deta_d};

    // Active-stress input of this quadrature point: BCS, or a prescribed
    // tension that replaces it (the BCS states then stay decoupled)
    ActiveInput ai;
    if (active_tension_.defined()) {
      ai.tension = tension;
    } else {
      ai.tauc = Y(i_tauc(q));
      ai.ecdot = Yd(i_ec(q));
      ai.mu = bp.mu;
    }

    // Hydrostatic pressure of this element (mixed u/p)
    p.p_mix = Y(i_pmix(pt.elem));

    // Internal forces g and their derivatives w.r.t. tau_c, p and the rates
    Kinematics<double> kin = compute_kinematics<double>(xi, pt.R);
    Matrix3d M = fiber_tensor<double>(pt.alpha);
    Matrix3d S = compute_stress<double>(kin, xidot, ai, pt.alpha, p);
    const double I4q = (kin.C * M).trace();  // e_F . C . e_F
    double g[6], dgdtau[6], dgdp[6], Kd[6][6], K[6][6];
    // d g_k / d p = (J C^-1) : dE_k
    const Matrix3d JCinv = (kin.J * kin.C.inverse()).eval();
    for (int i = 0; i < 6; i++) {
      g[i] = contract<double>(S, kin.dE[i]);
      dgdtau[i] = contract<double>(M, kin.dE[i]);
      dgdp[i] = contract<double>(JCinv, kin.dE[i]);
      for (int j = 0; j < 6; j++)
        Kd[i][j] = p.gamma * contract<double>(kin.dE[i], kin.dE[j]);
    }
    // Complex-step tangent K[k][m] = d g_k / d xi_m at fixed rates
    {
      const double h = 1e-30;
      cplx xic[6], xidotc[6], gc[6];
      for (int i = 0; i < 6; i++) {
        xic[i] = xi[i];
        xidotc[i] = xidot[i];
      }
      for (int m = 0; m < 6; m++) {
        xic[m] += cplx(0.0, h);
        compute_g<cplx>(xic, xidotc, ai, pt.R, pt.alpha, p, gc);
        for (int k = 0; k < 6; k++) K[k][m] = std::imag(gc[k]) / h;
        xic[m] = xi[m];
      }
    }

    // d sigma_1D / d tau_c = 1 / (1 + e_fib) (0 with a prescribed tension)
    const int active_col = i_tauc(q);
    const double i4denom = std::sqrt(I4q);  // 1 + e_fib = sqrt(I4)
    const double active_coef = active_tension_.defined() ? 0.0 : 1.0 / i4denom;

    // Equation rows and variable columns of each xi component with their
    // shape-function weights: 0 = rho (N), 1 = rho' (dN), 2 = phi' (dN),
    // 3 = beta, 4 = eps, 5 = eta' (dN).
    int req[6][2], cvr[6][2], rn[6], cn[6];
    double rw[6][2], cw[6][2];
    req[0][0] = e_rho(A0);
    rw[0][0] = N0;
    req[0][1] = e_rho(A1);
    rw[0][1] = N1;
    rn[0] = 2;
    cvr[0][0] = i_rho(A0);
    cw[0][0] = N0;
    cvr[0][1] = i_rho(A1);
    cw[0][1] = N1;
    cn[0] = 2;
    req[1][0] = e_rho(A0);
    rw[1][0] = dN0;
    req[1][1] = e_rho(A1);
    rw[1][1] = dN1;
    rn[1] = 2;
    cvr[1][0] = i_rho(A0);
    cw[1][0] = dN0;
    cvr[1][1] = i_rho(A1);
    cw[1][1] = dN1;
    cn[1] = 2;
    {  // phi': skip the row of the pinned inner node
      int c = 0;
      if (!pin_rigid_ || A0 != 0) {
        req[2][c] = e_phi(A0);
        rw[2][c] = dN0;
        c++;
      }
      req[2][c] = e_phi(A1);
      rw[2][c] = dN1;
      c++;
      rn[2] = c;
    }
    cvr[2][0] = i_phi(A0);
    cw[2][0] = dN0;
    cvr[2][1] = i_phi(A1);
    cw[2][1] = dN1;
    cn[2] = 2;
    req[3][0] = e_beta();
    rw[3][0] = 1.0;
    rn[3] = 1;
    cvr[3][0] = i_beta();
    cw[3][0] = 1.0;
    cn[3] = 1;
    req[4][0] = e_eps();
    rw[4][0] = 1.0;
    rn[4] = 1;
    cvr[4][0] = i_eps();
    cw[4][0] = 1.0;
    cn[4] = 1;
    {  // eta': skip the row of the pinned inner node
      int c = 0;
      if (!pin_rigid_ || A0 != 0) {
        req[5][c] = e_eta(A0);
        rw[5][c] = dN0;
        c++;
      }
      req[5][c] = e_eta(A1);
      rw[5][c] = dN1;
      c++;
      rn[5] = c;
    }
    cvr[5][0] = i_eta(A0);
    cw[5][0] = dN0;
    cvr[5][1] = i_eta(A1);
    cw[5][1] = dN1;
    cn[5] = 2;

    const double pref = 2.0 * M_PI * Lp * pt.w;
    for (int k = 0; k < 6; k++) {
      for (int ri = 0; ri < rn[k]; ri++) {
        const int re = req[k][ri];
        const double rwt = pref * rw[k][ri];
        Cloc[re] += rwt * g[k];
        Kat(re, active_col) += rwt * active_coef * dgdtau[k];
        // coupling to the element pressure p
        Kat(re, i_pmix(pt.elem)) += rwt * dgdp[k];
        // sigma_1D depends on e_c_dot through mu e_c_dot
        Kdat(re, i_ec(q)) += rwt * (ai.mu / i4denom) * dgdtau[k];
        for (int m = 0; m < 6; m++) {
          for (int ci = 0; ci < cn[m]; ci++) {
            const int cv = cvr[m][ci];
            Kat(re, cv) += rwt * K[k][m] * cw[m][ci];
            Kdat(re, cv) += rwt * Kd[k][m] * cw[m][ci];
          }
        }
      }
    }

    // --- Convective inertia D2u(zeta_dot, zeta_dot) (genet23 Eqs. 8, 10) ---
    // Acceleration a_c = 2 rho_dot w_om e_theta - r w_om^2 e_r with angular
    // velocity w_om = beta_dot Z + phi_dot and current radius r = R + rho,
    // tested with the virtual velocity rho^ e_r + (beta^ Z + phi^) r e_theta
    // and integrated over Theta and Z.
    {
      const double rho0 = parameters[global_param_ids[ParamId::density]];
      const double Wr = N0 * Y(i_vrho(A0)) + N1 * Y(i_vrho(A1));  // rho_dot
      const double Wp = N0 * Y(i_vphi(A0)) + N1 * Y(i_vphi(A1));  // phi_dot
      const double Wb = Y(i_vbeta());                             // beta_dot
      const double c0 = 2.0 * M_PI * rho0 * pt.w;
      const double base = c0 * (pt.R + rho);  // d(base)/d(rho) = c0
      const double L1 = Lp, L2 = Lp * Lp / 2.0, L3 = Lp * Lp * Lp / 3.0;
      const double Iw2 =
          Wb * Wb * L3 + 2.0 * Wb * Wp * L2 + Wp * Wp * L1;  // int w_om^2
      const double Qp = Wb * L2 + Wp * L1;                   // int w_om
      const double Qb = Wb * L3 + Wp * L2;                   // int w_om Z
      const int An[2] = {A0, A1};
      const double Nn[2] = {N0, N1};
      // beta (global twist): + 2 rho_dot (int w_om Z)
      Cloc[e_beta()] += base * 2.0 * Wr * Qb;
      Kat(e_beta(), i_vbeta()) += base * 2.0 * Wr * L3;
      for (int j = 0; j < 2; j++) {
        Kat(e_beta(), i_vrho(An[j])) += base * 2.0 * Nn[j] * Qb;
        Kat(e_beta(), i_vphi(An[j])) += base * 2.0 * Wr * L2 * Nn[j];
        Kat(e_beta(), i_rho(An[j])) += c0 * Nn[j] * 2.0 * Wr * Qb;
      }
      for (int i = 0; i < 2; i++) {
        const int a = An[i];
        const double Na = Nn[i];
        // rho: - (int w_om^2)
        Cloc[e_rho(a)] += -base * Iw2 * Na;
        Kat(e_rho(a), i_vbeta()) +=
            -base * Na * (2.0 * Wb * L3 + 2.0 * Wp * L2);
        for (int j = 0; j < 2; j++) {
          Kat(e_rho(a), i_vphi(An[j])) +=
              -base * Na * (2.0 * Wb * L2 + 2.0 * Wp * L1) * Nn[j];
          Kat(e_rho(a), i_rho(An[j])) += -c0 * Nn[j] * Iw2 * Na;
        }
        // phi (skip the inner node when pinned): + 2 rho_dot (int w_om)
        if (!pin_rigid_ || a != 0) {
          Cloc[e_phi(a)] += base * 2.0 * Wr * Qp * Na;
          Kat(e_phi(a), i_vbeta()) += base * 2.0 * Wr * L2 * Na;
          for (int j = 0; j < 2; j++) {
            Kat(e_phi(a), i_vrho(An[j])) += base * 2.0 * Nn[j] * Qp * Na;
            Kat(e_phi(a), i_vphi(An[j])) += base * 2.0 * Wr * L1 * Nn[j] * Na;
            Kat(e_phi(a), i_rho(An[j])) += c0 * Nn[j] * 2.0 * Wr * Qp * Na;
          }
        }
      }
    }

    // --- Current-radius correction of the mass (genet23 Eq. A5) ---
    // The beta and phi moment arms use the current radius: R^2 -> (R + rho)^2.
    // update_constant assembles the R^2 part; the difference is added to the
    // residual here. Its derivative w.r.t. rho is left out of the Jacobian
    // (inexact Newton, more robust); the residual and solution are exact.
    {
      const double rho0m = parameters[global_param_ids[ParamId::density]];
      const double cL1 = rho0m * 2.0 * M_PI * Lp;
      const double cL2 = rho0m * 2.0 * M_PI * Lp * Lp / 2.0;
      const double cL3 = rho0m * 2.0 * M_PI * Lp * Lp * Lp / 3.0;
      const double rcur = pt.R + rho;
      const double dr2 = rcur * rcur - pt.R * pt.R;  // (R+rho)^2 - R^2
      const int Am[2] = {A0, A1};
      const double Nm[2] = {N0, N1};
      // Accelerations beta_ddot and phi_ddot
      const double aBeta = Yd(i_vbeta());
      const double aPhi = N0 * Yd(i_vphi(A0)) + N1 * Yd(i_vphi(A1));
      // beta momentum row
      const double fBeta = cL3 * aBeta + cL2 * aPhi;
      Cloc[e_beta()] += pt.w * dr2 * fBeta;
      Kdat(e_beta(), i_vbeta()) += pt.w * dr2 * cL3;
      for (int j = 0; j < 2; j++)
        Kdat(e_beta(), i_vphi(Am[j])) += pt.w * dr2 * cL2 * Nm[j];
      // phi momentum rows (skip the pinned inner node)
      const double fPhi = cL2 * aBeta + cL1 * aPhi;
      for (int i = 0; i < 2; i++) {
        const int a = Am[i];
        if (pin_rigid_ && a == 0) continue;
        Cloc[e_phi(a)] += pt.w * dr2 * Nm[i] * fPhi;
        Kdat(e_phi(a), i_vbeta()) += pt.w * dr2 * Nm[i] * cL2;
        for (int j = 0; j < 2; j++)
          Kdat(e_phi(a), i_vphi(Am[j])) += pt.w * dr2 * Nm[i] * cL1 * Nm[j];
      }
    }

    // --- Mixed u/p: element pressure equation (Eq. 27) ---
    // Incompressible: int_e (J - 1) dOmega = 0. With a bulk modulus kappa
    // (ST91, p = kappa/2 (J - 1/J)): int_e ((J - 1/J)/2 - p/kappa) dOmega = 0.
    // J depends on rho, rho' and eps only.
    {
      const int erow = e_pmix(pt.elem);
      double dJfac = 1.0;  // d(constraint)/dJ
      if (kappa > 0.0) {
        Cloc[erow] += pref * (0.5 * (kin.J - 1.0 / kin.J) - p.p_mix / kappa);
        dJfac = 0.5 * (1.0 + 1.0 / (kin.J * kin.J));
        Kat(erow, i_pmix(pt.elem)) += -pref / kappa;
      } else {
        Cloc[erow] += pref * (kin.J - 1.0);
      }
      const double Rq = pt.R;
      double dJ[6] = {0, 0, 0, 0, 0, 0};
      dJ[0] = (1.0 + xi[1]) * (1.0 / Rq) * (1.0 + xi[4]);  // dJ/drho
      dJ[1] = (1.0 + xi[0] / Rq) * (1.0 + xi[4]);          // dJ/drho'
      dJ[4] = (1.0 + xi[1]) * (1.0 + xi[0] / Rq);          // dJ/deps
      for (int m = 0; m < 6; m++) {
        if (dJ[m] == 0.0) continue;
        for (int ci = 0; ci < cn[m]; ci++)
          Kat(erow, cvr[m][ci]) += pref * dJfac * dJ[m] * cw[m][ci];
      }
    }

    // --- BCS internal-variable ODEs at this quadrature point (Eqs. 32-33) ---
    {
      const double sfib = std::sqrt(I4q);
      const double e_fib = sfib - 1.0;
      const int erow[3] = {e_ec(q), e_tauc(q), e_kc(q)};
      const int avar[3] = {i_ec(q), i_tauc(q), i_kc(q)};
      double a[3] = {Y(i_ec(q)), Y(i_tauc(q)), Y(i_kc(q))};
      double ad[3] = {Yd(i_ec(q)), Yd(i_tauc(q)), Yd(i_kc(q))};
      double Ra[3];
      bcs_residual<double>(a, ad, e_fib, act, act_plus, bp, n0_ec_, n0_val_,
                           Ra);
      for (int i = 0; i < 3; i++) Cloc[erow[i]] += Ra[i];
      // Complex-step tangents w.r.t. the BCS states and their rates
      {
        const double h = 1e-30;
        cplx ac[3], adc[3], Rc[3];
        for (int t = 0; t < 3; t++) {
          ac[t] = a[t];
          adc[t] = ad[t];
        }
        for (int jj = 0; jj < 3; jj++) {
          ac[jj] += cplx(0.0, h);
          bcs_residual<cplx>(ac, adc, cplx(e_fib), act, act_plus, bp, n0_ec_,
                             n0_val_, Rc);
          for (int ii = 0; ii < 3; ii++)
            Kat(erow[ii], avar[jj]) += std::imag(Rc[ii]) / h;
          ac[jj] = a[jj];
          adc[jj] += cplx(0.0, h);
          bcs_residual<cplx>(ac, adc, cplx(e_fib), act, act_plus, bp, n0_ec_,
                             n0_val_, Rc);
          for (int ii = 0; ii < 3; ii++)
            Kdat(erow[ii], avar[jj]) += std::imag(Rc[ii]) / h;
          adc[jj] = ad[jj];
        }
      }
      // Only the e_c ODE depends on the wall, through e_fib:
      // d e_fib / d xi_m = dgdtau[m] / sqrt(I4).
      const double dRec_defib = -bp.k_s / bp.mu;
      for (int m = 0; m < 6; m++) {
        const double coef = dRec_defib * dgdtau[m] / sfib;
        for (int ci = 0; ci < cn[m]; ci++)
          Kat(e_ec(q), cvr[m][ci]) += coef * cw[m][ci];
      }
    }
  }

  // --- External pressure loading on the inner surface and lids (Eq. 63) ---
  const double rho_i = Y(i_rho(0));
  const double Rinner = Rip + rho_i;
  const double twoPiL = 2.0 * M_PI * Lp;
  const double piL = M_PI * Lp;
  // rho equation of the inner node: - P 2 pi L Rinner (1 + eps)
  {
    const int re = e_rho(0);
    Cloc[re] -= Pin * twoPiL * Rinner * (1.0 + eps);
    Kat(re, 0) -= twoPiL * Rinner * (1.0 + eps);      // d/dPin
    Kat(re, i_rho(0)) -= Pin * twoPiL * (1.0 + eps);  // d/drho_i
    Kat(re, i_eps()) -= Pin * twoPiL * Rinner;        // d/deps
  }
  // eps equation: - P pi L Rinner^2
  {
    const int re = e_eps();
    Cloc[re] -= Pin * piL * Rinner * Rinner;
    Kat(re, 0) -= piL * Rinner * Rinner;            // d/dPin
    Kat(re, i_rho(0)) -= Pin * piL * 2.0 * Rinner;  // d/drho_i
  }

  // --- Cavity volume V - pi Rinner^2 (1 + eps) L = 0 ---
  {
    const int re = e_volume();
    Cloc[re] -= piL * Rinner * Rinner * (1.0 + eps);
    Kat(re, i_rho(0)) -= piL * 2.0 * Rinner * (1.0 + eps);
    Kat(re, i_eps()) -= piL * Rinner * Rinner;
  }

  // --- Mass conservation Qin - Qout - Vdot = 0 ---
  {
    const int re = e_mass();
    const double rho_i_d = Yd(i_rho(0));
    const double Vdot =
        piL * (2.0 * Rinner * (1.0 + eps) * rho_i_d + Rinner * Rinner * eps_d);
    Cloc[re] -= Vdot;
    // d/dy
    Kat(re, i_rho(0)) -=
        piL * (2.0 * (1.0 + eps) * rho_i_d + 2.0 * Rinner * eps_d);
    Kat(re, i_eps()) -= piL * 2.0 * Rinner * rho_i_d;
    // d/dydot
    Kdat(re, i_rho(0)) -= piL * 2.0 * Rinner * (1.0 + eps);
    Kdat(re, i_eps()) -= piL * Rinner * Rinner;
  }

  // Fixed length: the eps row holds only the constraint eps = 0 (in F).
  if (fixed_length_) {
    const int re = e_eps();
    Cloc[re] = 0.0;
    for (int j = 0; j < n_var; j++) Kat(re, j) = Kdat(re, j) = 0.0;
  }

  // --- Assign the local contributions to the sparse system ---
  for (int i = 0; i < n_eqn; i++) {
    system.C(global_eqn_ids[i]) = Cloc[i];
    for (int j = 0; j < n_var; j++) {
      system.dC_dy.coeffRef(global_eqn_ids[i], global_var_ids[j]) = Kat(i, j);
      system.dC_dydot.coeffRef(global_eqn_ids[i], global_var_ids[j]) =
          Kdat(i, j);
    }
  }
}

void ChamberCylinder::set_frank_starling_curve(std::vector<double> ec,
                                               std::vector<double> val) {
  if (ec.size() != val.size() || ec.size() < 2) {
    throw std::runtime_error(
        "ChamberCylinder::set_frank_starling_curve: e_c and n0 arrays must "
        "have equal length >= 2 (got " +
        std::to_string(ec.size()) + " and " + std::to_string(val.size()) +
        ").");
  }
  n0_ec_ = std::move(ec);
  n0_val_ = std::move(val);
}

void FourierTension::set(const std::vector<double>& t,
                         const std::vector<double>& v, int n_modes,
                         bool use_ramp) {
  const size_t n = t.size();
  if (n < 2 || v.size() != n || n_modes < 1) {
    throw std::runtime_error(
        "FourierTension: needs >= 2 samples t, T of equal length and "
        "fourier_modes >= 1.");
  }
  for (size_t i = 1; i < n; i++) {
    if (t[i] <= t[i - 1]) {
      throw std::runtime_error(
          "FourierTension: sample times must increase strictly.");
    }
  }
  ramp = use_ramp;
  t0 = t[0];
  period = t[n - 1] - t[0];
  v0 = v[0];
  slope = (v[n - 1] - v[0]) / period;
  // Samples minus the linear trend
  std::vector<double> ts(n), vs(n);
  for (size_t i = 0; i < n; i++) {
    ts[i] = t[i] - t0;
    vs[i] = v[i] - (v0 + slope * ts[i]);
  }
  // Exact Fourier coefficients of the piecewise-linear interpolant
  const int n_coef = ramp ? 1 : n_modes;
  re.assign(n_coef, 0.0);
  im.assign(n_coef, 0.0);
  for (int k = 0; k < n_coef; k++) {
    for (size_t i = 0; i + 1 < n; i++) {
      if (k == 0) {
        re[k] += 0.5 * (ts[i + 1] - ts[i]) * (vs[i + 1] + vs[i]);
      } else {
        const double ko = 2.0 * M_PI * k * ts[i] / period;
        const double kn = 2.0 * M_PI * k * ts[i + 1] / period;
        const double s = (vs[i + 1] - vs[i]) / (ts[i + 1] - ts[i]);
        re[k] += s * (std::cos(kn) - std::cos(ko));
        im[k] -= s * (std::sin(kn) - std::sin(ko));
      }
    }
    if (k == 0) {
      re[k] /= period;
    } else {
      const double scale = 0.5 * period / (M_PI * M_PI * k * k);
      re[k] *= scale;
      im[k] *= scale;
    }
  }
}

double FourierTension::value(double time) const {
  const double t = ramp ? std::max(std::min(time - t0, period), 0.0)
                        : std::fmod(time - t0, period);
  double val = v0 + t * slope;
  if (!ramp) {
    for (size_t k = 0; k < re.size(); k++) {
      const double K = 2.0 * M_PI * k * t / period;
      val += re[k] * std::cos(K) - im[k] * std::sin(K);
    }
  }
  return val;
}
