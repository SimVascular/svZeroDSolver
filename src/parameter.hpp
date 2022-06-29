#ifndef SVZERODSOLVER_PARAMETER_H_
#define SVZERODSOLVER_PARAMETER_H_

#include "dofhandler.hpp"
#include <math.h>

#include <vector>
#include <numeric>

template <typename T>
class TimeDependentParameter
{
public:
    TimeDependentParameter();
    TimeDependentParameter(std::vector<T> times, std::vector<T> values);
    ~TimeDependentParameter();

    std::vector<T> times;
    std::vector<T> values;
    T cycle_period;
    int size;

    T get(T time);

    void setup_dofs(DOFHandler &dofhandler);
    bool isconstant;

    void to_steady();
};

template <typename T>
TimeDependentParameter<T>::TimeDependentParameter()
{
}

template <typename T>
TimeDependentParameter<T>::TimeDependentParameter(std::vector<T> times, std::vector<T> values)
{
    this->times = times;
    this->values = values;
    cycle_period = times.back() - times[0];
    size = times.size();
    if (size == 1)
    {
        isconstant = true;
    }
    else
    {
        isconstant = false;
    }
}

template <typename T>
TimeDependentParameter<T>::~TimeDependentParameter()
{
}

template <typename T>
T TimeDependentParameter<T>::get(T time)
{
    if (size == 1)
    {
        return values[0];
    }

    T rtime = fmod(time, cycle_period);

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

template <typename T>
void TimeDependentParameter<T>::to_steady()
{
    T mean = std::accumulate(values.begin(), values.end(), 0.0) / T(values.size());
    values = std::vector<T>();
    times = std::vector<T>();
    values.push_back(mean);
    size = 1;
    isconstant = true;
}

#endif // SVZERODSOLVER_PARAMETER_H_