
#include"PressureReferenceBC.h"

namespace zd_model {

void PressureReferenceBC::setup_dofs(DOFHandler &dofhandler) 
{
  Block::setup_dofs_(dofhandler, 1, {});
}

void PressureReferenceBC::update_constant(algebra::SparseSystem &system,
                                             std::vector<double> &parameters) 
{
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[0]) = 1.0;
}

void PressureReferenceBC::update_time(algebra::SparseSystem &system,
                                         std::vector<double> &parameters) 
{
  system.C(this->global_eqn_ids[0]) = -parameters[this->global_param_ids[0]];
}

std::map<std::string, int> PressureReferenceBC::get_num_triplets() {
  return num_triplets;
}

}  
