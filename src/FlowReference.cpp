#include "flowreference.hpp"

void FlowReference::update_constant(System system)
{
    system.F(global_eqn_ids[0], global_var_ids[1]) = 1.0;
}
void FlowReference::update_time(System system)
{
    system.C(global_eqn_ids[0]) = -params->Q;
}