#ifndef SVZERODSOLVER_MODEL_H_
#define SVZERODSOLVER_MODEL_H_

#include <map>
#include <list>

#include "junction.hpp"
#include "bloodvessel.hpp"
#include "rcrblockwithdistalpressure.hpp"
#include "flowreference.hpp"
#include "dofhandler.hpp"
#include "node.hpp"

class Model
{
public:
    Model();
    ~Model();

    std::map<std::string, std::variant<Junction, BloodVessel, FlowReference, RCRBlockWithDistalPressure>> blocks;
    DOFHandler dofhandler;
    std::list<Node> nodes;
};

#endif // SVZERODSOLVER_MODEL_H_