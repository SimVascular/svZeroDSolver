// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#ifndef SVZERODSOLVER_MODEL_CLOSEDLOOPHEARTPULMONARYSMOOTH_HPP_
#define SVZERODSOLVER_MODEL_CLOSEDLOOPHEARTPULMONARYSMOOTH_HPP_

#include "Block.h"
#include "SparseSystem.h"

#define _USE_MATH_DEFINES
#include <cmath>

/**
 * @brief Monolithic heart+pulmonary with smooth (tanh) valves.
 *
 * Identical to ClosedLoopHeartPulmonary except the four valve positions
 * use a smooth sigmoid s(ΔP) = 0.5*(1 + tanh(k*ΔP)) instead of hard
 * 0/1 switching. No post_solve Q zeroing. All other equations, activation
 * functions, and P-V relationships are unchanged.
 *
 * Additional parameter: Steepness (valve sigmoid sharpness).
 */
class ClosedLoopHeartPulmonarySmooth : public Block {
 public:
  ClosedLoopHeartPulmonarySmooth(int id, Model* model)
      : Block(id, model, BlockType::closed_loop_heart_pulmonary_smooth,
              BlockClass::closed_loop,
              {{"Tsa", InputParameter()},     {"tpwave", InputParameter()},
               {"Erv_s", InputParameter()},   {"Elv_s", InputParameter()},
               {"iml", InputParameter()},     {"imr", InputParameter()},
               {"Lra_v", InputParameter()},   {"Rra_v", InputParameter()},
               {"Lrv_a", InputParameter()},   {"Rrv_a", InputParameter()},
               {"Lla_v", InputParameter()},   {"Rla_v", InputParameter()},
               {"Llv_a", InputParameter()},   {"Rlv_ao", InputParameter()},
               {"Vrv_u", InputParameter()},   {"Vlv_u", InputParameter()},
               {"Rpd", InputParameter()},     {"Cp", InputParameter()},
               {"Cpa", InputParameter()},     {"Kxp_ra", InputParameter()},
               {"Kxv_ra", InputParameter()},  {"Kxp_la", InputParameter()},
               {"Kxv_la", InputParameter()},  {"Emax_ra", InputParameter()},
               {"Emax_la", InputParameter()}, {"Vaso_ra", InputParameter()},
               {"Vaso_la", InputParameter()}, {"Steepness", InputParameter()},
               {"Rmax", InputParameter()}}) {}

  enum ParamId {
    TSA = 0, TPWAVE = 1, ERV_S = 2, ELV_S = 3, IML = 4, IMR = 5,
    LRA_V = 6, RRA_V = 7, LRV_A = 8, RRV_A = 9, LLA_V = 10, RLA_V = 11,
    LLV_A = 12, RLV_AO = 13, VRV_U = 14, VLV_U = 15, RPD = 16, CP = 17,
    CPA = 18, KXP_RA = 19, KXV_RA = 20, KXP_LA = 21, KXV_LA = 22,
    EMAX_RA = 23, EMAX_LA = 24, VASO_RA = 25, VASO_LA = 26,
    STEEPNESS = 27,
    RMAX = 28
  };

  void setup_dofs(DOFHandler& dofhandler);
  void update_constant(SparseSystem& system, std::vector<double>& parameters);
  void update_time(SparseSystem& system, std::vector<double>& parameters);
  void update_solution(SparseSystem& system, std::vector<double>& parameters,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
                       const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy);

  TripletsContributions num_triplets{33, 10, 2};

 private:
  double AA = 0.0;
  double Elv = 0.0;
  double Erv = 0.0;
  double psi_ra, psi_la, psi_ra_derivative, psi_la_derivative;
  double valves[16];       // sigmoid gate [0,1] for volume equations
  double r_scale[16];      // resistance scale for outflow equations
  double last_time_ = -1;

  void get_activation_and_elastance_functions(std::vector<double>& parameters);
  void get_psi_ra_la(std::vector<double>& parameters,
                     const Eigen::Matrix<double, Eigen::Dynamic, 1>& y);
  void get_valve_positions(std::vector<double>& parameters,
                           const Eigen::Matrix<double, Eigen::Dynamic, 1>& y);
};

#endif
