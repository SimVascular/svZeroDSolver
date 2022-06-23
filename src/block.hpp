#ifndef SVZERODSOLVER_BLOCK_H_
#define SVZERODSOLVER_BLOCK_H_

#include <vector>
#include "node.hpp"
#include "dofhandler.hpp"
#include "system.hpp"

class Block
{
protected:
    unsigned short int num_equations;
    unsigned short int num_internal_vars;

public:
    struct Parameters
    {
    };
    Block(Parameters &params, std::string name);
    ~Block();

    std::string name;
    std::vector<Node *> inlet_nodes;
    std::vector<Node *> outlet_nodes;

    std::vector<unsigned int> global_var_ids;
    std::vector<unsigned int> global_eqn_ids;

    void setup_dofs(DOFHandler &dofhandler);
    void update_constant(System system);
    void update_time(System system, double time);
    void update_solution(System system, Eigen::VectorXd &y);

private:
    Parameters *params;
};

#endif // SVZERODSOLVER_BLOCK_H_