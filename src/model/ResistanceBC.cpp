
#include "ResistanceBC.h"

namespace zd_model {

void ResistanceBC::setup_dofs(DOFHandler &dofhandler) 
{
  Block::setup_dofs_(dofhandler, 1, {});
}

void ResistanceBC::update_constant(algebra::SparseSystem &system, std::vector<double> &parameters) 
{
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[0]) = 1.0;
}

void ResistanceBC::update_time(algebra::SparseSystem &system, std::vector<double> &parameters) 
{
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) =
      -parameters[this->global_param_ids[0]];
  system.C(this->global_eqn_ids[0]) = -parameters[this->global_param_ids[1]];
}

std::map<std::string, int> ResistanceBC::get_num_triplets() {
  return num_triplets;
}

}  
