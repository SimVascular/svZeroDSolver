#ifndef SVZERODSOLVER_MODEL_H_
#define SVZERODSOLVER_MODEL_H_

#include <map>

#include "junction.hpp"
#include "bloodvessel.hpp"
#include "rcrblockwithdistalpressure.hpp"
#include "flowreference.hpp"
#include "dofhandler.hpp"

class Model
{
public:
    Model();
    ~Model();

    std::map<std::string, std::variant<Junction, BloodVessel, FlowReference, RCRBlockWithDistalPressure>> blocks;
    DOFHandler dofhandler;
};

#endif // SVZERODSOLVER_MODEL_H_