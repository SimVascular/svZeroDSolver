
#include "ClosedLoopCoronaryBC.h"
#include "Model.h"

namespace zd_model {

void ClosedLoopCoronaryBC::setup_dofs(DOFHandler &dofhandler) {
  Block::setup_dofs_(dofhandler, 3, {"volume_im"});
}

void ClosedLoopCoronaryBC::setup_model_dependent_params() 
{
  auto heart_block = this->model->get_block("CLH");

  if (side == Side::LEFT) {
    im_param_id = heart_block->global_param_ids[ClosedLoopHeartPulmonary::ParamId::IML];
    ventricle_var_id = heart_block->global_var_ids[13];  // Solution ID for LV pressure

  } else {
    im_param_id = heart_block->global_param_ids[ClosedLoopHeartPulmonary::ParamId::IMR];
    ventricle_var_id = heart_block->global_var_ids[6];
  }
}

void ClosedLoopCoronaryBC::update_constant(algebra::SparseSystem &system, std::vector<double> &parameters) 
{
  auto ra = parameters[this->global_param_ids[ParamId::RA]];
  auto ram = parameters[this->global_param_ids[ParamId::RAM]];
  auto rv = parameters[this->global_param_ids[ParamId::RV]];
  auto ca = parameters[this->global_param_ids[ParamId::CA]];
  auto cim = parameters[this->global_param_ids[ParamId::CIM]];

  system.E.coeffRef(this->global_eqn_ids[0], this->global_var_ids[0]) =
      -ram * ca;
  system.E.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) =
      ram * ra * ca;
  system.E.coeffRef(this->global_eqn_ids[1], this->global_var_ids[0]) = -ca;
  system.E.coeffRef(this->global_eqn_ids[1], this->global_var_ids[1]) = ca * ra;
  system.E.coeffRef(this->global_eqn_ids[1], this->global_var_ids[4]) = -1.0;

  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[0]) = -1.0;
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) =
      (ra + ram);
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[2]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[3]) = rv;
  system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[1]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[3]) = -1.0;
  system.F.coeffRef(this->global_eqn_ids[2], this->global_var_ids[2]) = cim;
  system.F.coeffRef(this->global_eqn_ids[2], this->global_var_ids[3]) =
      cim * rv;
  system.F.coeffRef(this->global_eqn_ids[2], this->global_var_ids[4]) = -1.0;
}

void ClosedLoopCoronaryBC::update_solution(
    algebra::SparseSystem &system, std::vector<double> &parameters,
    Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
    Eigen::Matrix<double, Eigen::Dynamic, 1> &dy) 
{
  auto cim = parameters[this->global_param_ids[ParamId::CIM]];
  auto im = parameters[im_param_id];
  auto pim = im * y[this->ventricle_var_id];
  system.C(this->global_eqn_ids[2]) = -cim * pim;
}

std::map<std::string, int> ClosedLoopCoronaryBC::get_num_triplets() 
{
  return num_triplets;
}

}  
