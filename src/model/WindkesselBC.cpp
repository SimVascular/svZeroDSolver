
#include "WindkesselBC.h"

namespace zd_model {

void WindkesselBC::setup_dofs(DOFHandler &dofhandler) {
  Block::setup_dofs_(dofhandler, 2, {"pressure_c"});
}

void WindkesselBC::update_constant(algebra::SparseSystem &system,
                                      
std::vector<double> &parameters) {
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[0]) = 1.0;
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[2]) = -1.0;
  system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[2]) = -1.0;
}

void WindkesselBC::update_time(algebra::SparseSystem &system,
                                  
std::vector<double> &parameters) {
  system.E.coeffRef(this->global_eqn_ids[1], this->global_var_ids[2]) =
      -parameters[this->global_param_ids[2]] *
      parameters[this->global_param_ids[1]];
  system.F.coeffRef(this->global_eqn_ids[0], this->global_var_ids[1]) =
      -parameters[this->global_param_ids[0]];
  system.F.coeffRef(this->global_eqn_ids[1], this->global_var_ids[1]) =
      parameters[this->global_param_ids[2]];
  system.C(this->global_eqn_ids[1]) = parameters[this->global_param_ids[3]];
}

std::map<std::string, int> WindkesselBC::get_num_triplets() {
  return num_triplets;
}

}  
