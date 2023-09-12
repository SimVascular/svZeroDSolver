
#include "ClosedLoopRCRBC.h"

namespace zd_model {

void ClosedLoopRCRBC::setup_dofs(DOFHandler &dofhandler) 
{
  Block::setup_dofs_(dofhandler, 3, {"P_c"});
}

void ClosedLoopRCRBC::update_constant(algebra::SparseSystem &system, std::vector<double> &parameters) 
{
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) = -1.0;
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[3]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[0]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[4]) = -1.0;
  system.F.coeffRef(this->global_eqn_ids[2], this->global_var_ids[2]) = -1.0;
  system.F.coeffRef(this->global_eqn_ids[2], this->global_var_ids[4]) = 1.0;

  // Below values can be unsteady if needed (not currently implemented)
  system.E.coeffRef(this->global_eqn_ids[0], this->global_var_ids[4]) =
      parameters[this->global_param_ids[ParamId::C]];
  system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[1]) =
      -parameters[this->global_param_ids[ParamId::RP]];
  system.F.coeffRef(this->global_eqn_ids[2], this->global_var_ids[3]) =
      -parameters[this->global_param_ids[ParamId::RD]];
}

std::map<std::string, int> ClosedLoopRCRBC::get_num_triplets() 
{
  return num_triplets;
}

}  
