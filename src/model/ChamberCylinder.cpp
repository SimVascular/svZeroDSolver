// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ChamberCylinder.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <initializer_list>
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

/// Material parameters passed to the continuum kernel.
struct MatParams {
  const Material* material = nullptr;  ///< wall material (deviatoric stress)
  double gamma = 0.0;                  ///< viscosity
  double p_mix = 0.0;  ///< mixed u/p element hydrostatic pressure
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

/**
 * @brief Second Piola-Kirchhoff stress at an integration point.
 *
 * Sigma = Sigma^d (material, \cite genet23 Eq. 25) + Sigma^b (hydrostatic
 * pressure, Eqs. 26-27) + Sigma^v (viscous, Eq. 28) + Sigma^a (active stress
 * law, Eqs. 29-30).
 */
template <typename T>
Mat3<T> compute_stress(const Kinematics<T>& k, const T xidot[6], double tau,
                       double alpha, const MatParams& p) {
  const Mat3<T>& C = k.C;
  const T J = k.J;
  const Mat3<T> M = fiber_tensor<T>(alpha);
  // fiber (0, cos alpha, sin alpha) in the (e_R, e_Th, e_Z) basis
  const Material::Vec3 f(0.0, std::cos(alpha), std::sin(alpha));

  const Mat3<T> Sd = p.material->stress(C, J, f);  // deviatoric

  // Hydrostatic stress p J C^-1 with the element pressure p (mixed u/p).
  const T bulk = T(p.p_mix);
  Mat3<T> Sb = bulk * J * C.inverse();  // bulk

  Mat3<T> Edot = Mat3<T>::Zero();  // viscous
  for (int m = 0; m < 6; m++) Edot += k.dE[m] * xidot[m];
  Mat3<T> Sv = p.gamma * Edot;

  Mat3<T> Sa = tau * M;  // active, tau of the active stress law

  return Sd + Sb + Sv + Sa;
}

/// Generalized internal forces g_k = Sigma : dE/dxi_k at an integration point.
template <typename T>
void compute_g(const T xi[6], const T xidot[6], double tau, double R,
               double alpha, const MatParams& p, T g[6]) {
  Kinematics<T> k = compute_kinematics<T>(xi, R);
  Mat3<T> S = compute_stress<T>(k, xidot, tau, alpha, p);
  for (int i = 0; i < 6; i++) g[i] = contract<T>(S, k.dE[i]);
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
  check_inputs();
  n_law = active_stress_->num_states();

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
  // Active stress states per quadrature point (state 0: tau)
  const std::vector<std::string> law_states = active_stress_->state_names();
  for (int q = 0; q < n_quad; q++) {
    for (const auto& name : law_states)
      int_names.push_back(name + "_" + std::to_string(q));
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

  // Local tangents; their sparsity pattern is recorded by the first assembly
  K_.assign(static_cast<size_t>(n_eqn) * n_var, 0.0);
  Kd_.assign(K_.size(), 0.0);
  in_K_.assign(K_.size(), 0);
  in_Kd_.assign(K_.size(), 0);
  pattern_K_.clear();
  pattern_Kd_.clear();
  pattern_recorded_ = false;

  // Reservation estimates: per element its kinematic, velocity, pressure and
  // active stress unknowns, plus the global rows and columns (beta, eps,
  // volume, pressure); few F/E entries
  const int n_local = 6 + 6 + 4 + 1 + 3 * n_law;
  num_triplets.D = n_ele * n_local * n_local + 8 * (n_eqn + n_var);
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
    auto Madd = [&](int r, int c, double v) { M[{r, c}] += v; };
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
  (void)system;
  (void)parameters;
  active_stress_->update_time(activation_function_->compute(model->time));
}

void ChamberCylinder::update_solution(
    SparseSystem& system, std::vector<double>& parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) {
  MatParams p;
  p.material = material_.get();
  p.gamma = parameters[global_param_ids[ParamId::gamma]];
  const double Rip = parameters[global_param_ids[ParamId::Ri]];
  const double Lp = parameters[global_param_ids[ParamId::length]];

  // Block DOFs and their rates
  auto Y = [&](int lv) { return y[global_var_ids[lv]]; };
  auto Yd = [&](int lv) { return dy[global_var_ids[lv]]; };
  const double beta = Y(i_beta()), eps = Y(i_eps());
  const double beta_d = Yd(i_beta()), eps_d = Yd(i_eps());
  const double Pin = Y(0);

  // Local residual and tangents, scattered into the system at the end. The
  // first assembly records the touched tangent entries (the sparsity
  // pattern); later assemblies must stay within it.
  std::vector<double> Cloc(n_eqn, 0.0);
  const bool record = !pattern_recorded_;
  auto touch = [&](int r, int c, std::vector<char>& in,
                   std::vector<int>& pattern) {
    const int i = r * n_var + c;
    if (!in[i]) {
      if (!record) {
        throw std::runtime_error(
            "ChamberCylinder block '" + get_name() +
            "': tangent entry outside the recorded sparsity pattern.");
      }
      in[i] = 1;
      pattern.push_back(i);
    }
    return i;
  };
  auto Kat = [&](int r, int c) -> double& {
    return K_[touch(r, c, in_K_, pattern_K_)];
  };
  auto Kdat = [&](int r, int c) -> double& {
    return Kd_[touch(r, c, in_Kd_, pattern_Kd_)];
  };

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

    // Active fiber stress tau of this quadrature point (state 0 of the law)
    const double tau = Y(i_law(q, 0));

    // Hydrostatic pressure of this element (mixed u/p)
    p.p_mix = Y(i_pmix(pt.elem));

    // Internal forces g and their derivatives w.r.t. tau_c, p and the rates
    Kinematics<double> kin = compute_kinematics<double>(xi, pt.R);
    Matrix3d M = fiber_tensor<double>(pt.alpha);
    Matrix3d S = compute_stress<double>(kin, xidot, tau, pt.alpha, p);
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
        compute_g<cplx>(xic, xidotc, tau, pt.R, pt.alpha, p, gc);
        for (int k = 0; k < 6; k++) K[k][m] = std::imag(gc[k]) / h;
        xic[m] = xi[m];
      }
    }

    // The wall depends on the active stress law through tau (coefficient 1)
    const int active_col = i_law(q, 0);
    const double lam = std::sqrt(I4q);  // fiber stretch

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
        Kat(re, active_col) += rwt * dgdtau[k];
        // coupling to the element pressure p
        Kat(re, i_pmix(pt.elem)) += rwt * dgdp[k];
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

    // --- Mixed u/p: element incompressibility int_e (J - 1) dOmega = 0 ---
    // (genet23 Eq. 27); J depends on rho, rho' and eps only.
    {
      const int erow = e_pmix(pt.elem);
      Cloc[erow] += pref * (kin.J - 1.0);
      const double Rq = pt.R;
      double dJ[6] = {0, 0, 0, 0, 0, 0};
      dJ[0] = (1.0 + xi[1]) * (1.0 / Rq) * (1.0 + xi[4]);  // dJ/drho
      dJ[1] = (1.0 + xi[0] / Rq) * (1.0 + xi[4]);          // dJ/drho'
      dJ[4] = (1.0 + xi[1]) * (1.0 + xi[0] / Rq);          // dJ/deps
      for (int m : {0, 1, 4}) {  // dJ/drho, dJ/drho', dJ/deps
        for (int ci = 0; ci < cn[m]; ci++)
          Kat(erow, cvr[m][ci]) += pref * dJ[m] * cw[m][ci];
      }
    }

    // --- Active stress law at this quadrature point ---
    {
      std::vector<double> a(n_law), ad(n_law);
      for (int i = 0; i < n_law; i++) {
        a[i] = Y(i_law(q, i));
        ad[i] = Yd(i_law(q, i));
      }
      ActiveStress::Linearization lin;
      active_stress_->linearize(a.data(), ad.data(), lam, lin);
      for (int i = 0; i < n_law; i++) {
        const int erow = e_law(q, i);
        Cloc[erow] += lin.R[i];
        for (int j = 0; j < n_law; j++) {
          Kat(erow, i_law(q, j)) += lin.dR_da[i * n_law + j];
          Kdat(erow, i_law(q, j)) += lin.dR_dadot[i * n_law + j];
        }
        // the law depends on the wall through the fiber stretch:
        // d lambda / d xi_m = (e_F . dE_m . e_F) / lambda = dgdtau[m] / lambda
        for (int m = 0; m < 6; m++) {
          const double coef = lin.dR_dlambda[i] * dgdtau[m] / lam;
          for (int ci = 0; ci < cn[m]; ci++)
            Kat(erow, cvr[m][ci]) += coef * cw[m][ci];
        }
      }
    }
  }

  // --- External pressure loading on the inner surface (Eq. 63) ---
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
  // --- Lids loaded by the cavity pressure; eps equation: - P pi L Rinner^2 ---
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

  // --- Assign the local contributions to the sparse system ---
  for (int i = 0; i < n_eqn; i++) system.C(global_eqn_ids[i]) = Cloc[i];
  if (record) {
    // Insert the pattern column by column (cheap insertion into the
    // column-major sparse matrices)
    auto by_column = [&](int a, int b) {
      const int ca = global_var_ids[a % n_var], cb = global_var_ids[b % n_var];
      if (ca != cb) return ca < cb;
      return global_eqn_ids[a / n_var] < global_eqn_ids[b / n_var];
    };
    std::sort(pattern_K_.begin(), pattern_K_.end(), by_column);
    std::sort(pattern_Kd_.begin(), pattern_Kd_.end(), by_column);
    pattern_recorded_ = true;
  }
  // Write the pattern entries (zero where a coupling vanishes at this state)
  // and reset the local arrays
  for (const int i : pattern_K_) {
    system.dC_dy.coeffRef(global_eqn_ids[i / n_var],
                          global_var_ids[i % n_var]) = K_[i];
    K_[i] = 0.0;
  }
  for (const int i : pattern_Kd_) {
    system.dC_dydot.coeffRef(global_eqn_ids[i / n_var],
                             global_var_ids[i % n_var]) = Kd_[i];
    Kd_[i] = 0.0;
  }
}

void ChamberCylinder::check_inputs() const {
  auto fail = [&](const std::string& msg) {
    throw std::runtime_error("ChamberCylinder '" + get_name() + "': " + msg);
  };
  if (!material_) fail("material not set. Provide a \"material\" JSON block.");
  if (!active_stress_)
    fail("active_stress not set. Provide an \"active_stress\" JSON block.");
  if (!activation_function_)
    fail(
        "activation_function not set. Provide an \"activation_function\" "
        "JSON block.");
  if (activation_function_->signal_type() !=
      active_stress_->required_signal_type())
    fail(std::string("active_stress law requires a ") +
         (active_stress_->required_signal_type() == ActivationSignalType::rate
              ? "rate-type (e.g. piecewise_rate)"
              : "normalized (e.g. half_cosine)") +
         " activation function.");
}
