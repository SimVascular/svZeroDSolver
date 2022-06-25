#ifndef SVZERODSOLVER_PARAMETER_H_
#define SVZERODSOLVER_PARAMETER_H_

#include "dofhandler.hpp"
#include "json.h"
#include <math.h>

#include <vector>

class TimeDependentParameter
{
public:
    TimeDependentParameter();
    TimeDependentParameter(std::vector<double> times, std::vector<double> values);
    ~TimeDependentParameter();

    std::vector<double> times;
    std::vector<double> values;
    double cycle_period;
    int size;

    double get(double time);

    void setup_dofs(DOFHandler &dofhandler);
    bool isconstant;
};

#endif // SVZERODSOLVER_PARAMETER_H_