#include "parameter.hpp"

#include <iostream>

TimeDependentParameter::TimeDependentParameter()
{
}

TimeDependentParameter::TimeDependentParameter(std::vector<double> times, std::vector<double> values)
{
    this->times = times;
    this->values = values;
    cycle_period = times.back() - times[0];
    size = times.size();
}

TimeDependentParameter::~TimeDependentParameter()
{
}

double TimeDependentParameter::get(double time)
{
    if (size == 1)
    {
        return values[0];
    }

    double rtime = remainder(time, cycle_period);

    auto i = lower_bound(times.begin(), times.end(), rtime);
    unsigned int k = i - times.begin();

    if (i == times.end())
        --i;
    else if (*i == rtime)
    {
        return values[k];
    }

    unsigned int l = k ? k - 1 : 1;

    return values[l] + ((values[k] - values[l]) / (times[k] - times[l])) * (rtime - times[l]);
}