
#include "OpenLoopCoronaryBC.h"

namespace zd_model {

void OpenLoopCoronaryBC::setup_dofs(DOFHandler &dofhandler) 
{
  Block::setup_dofs_(dofhandler, 2, {"volume_im"});
}

void OpenLoopCoronaryBC::update_constant(algebra::SparseSystem &system, std::vector<double> &parameters) 
{
  auto Ra = parameters[this->global_param_ids[0]];
  auto Ram = parameters[this->global_param_ids[1]];
  auto Rv = parameters[this->global_param_ids[2]];
  auto Ca = parameters[this->global_param_ids[3]];
  auto Cim = parameters[this->global_param_ids[4]];

  if (this->steady) {
    // Different assmembly for steady block to avoid singular system
    // and solve for the internal variable V_im inherently
    system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[0]) = -Cim;
    system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) =
        Cim * (Ra + Ram);
    system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[2]) = 1.0;
    system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[0]) = -1.0;
    system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[1]) =
        Ra + Ram + Rv;
  } else {
    system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) =
        Cim * Rv;
    system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[2]) = -1.0;
    system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[0]) =
        Cim * Rv;
    system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[1]) =
        -Cim * Rv * Ra;
    system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[2]) =
        -(Rv + Ram);

    system.E.coeffRef(this->global_eqn_ids[0], this->global_var_ids[0]) =
        -Ca * Cim * Rv;
    system.E.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) =
        Ra * Ca * Cim * Rv;
    system.E.coeffRef(this->global_eqn_ids[0], this->global_var_ids[2]) =
        -Cim * Rv;
    system.E.coeffRef(this->global_eqn_ids[1], this->global_var_ids[2]) =
        -Cim * Rv * Ram;
  }
}

void OpenLoopCoronaryBC::update_time(algebra::SparseSystem &system, std::vector<double> &parameters) 
{
  auto Ram = parameters[this->global_param_ids[1]];
  auto Rv = parameters[this->global_param_ids[2]];
  auto Cim = parameters[this->global_param_ids[4]];
  auto Pim = parameters[this->global_param_ids[5]];
  auto Pv = parameters[this->global_param_ids[6]];

  if (this->steady) {
    system.C(this->global_eqn_ids[0]) = -Cim * Pim;
    system.C(this->global_eqn_ids[1]) = Pv;
  } else {
    system.C(this->global_eqn_ids[0]) = Cim * (-Pim + Pv);
    system.C(this->global_eqn_ids[1]) =
        -Cim * (Rv + Ram) * Pim + Ram * Cim * Pv;
  }
}

std::map<std::string, int> OpenLoopCoronaryBC::get_num_triplets() 
{
  return num_triplets;
}

} 
