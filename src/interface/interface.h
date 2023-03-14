
#include "algebra/integrator.hpp"
#include "algebra/sparsesystem.hpp"
#include "algebra/state.hpp"
#include "helpers/debug.hpp"
#include "helpers/endswith.hpp"
#include "io/configreader.hpp"
#include "io/csvwriter.hpp"
#include "model/model.hpp"

#include <map>
#include <string>
#include <vector>

typedef double T;
template <typename TT>
using S = ALGEBRA::SparseSystem<TT>;

/**
 * @brief Interface class for calling svZeroD from external programs
 */

class SolverInterface 
{
  public:
    SolverInterface(const std::string& input_file_name);
    ~SolverInterface(); 

    static int problem_id_count_;
    static std::map<int,SolverInterface*> interface_list_; 

    int problem_id_ = 0;
    std::string input_file_name_;

    // Parameters for the external solver (the calling program).
    // This is set by the external solver via the interface.
    double external_step_size_ = 0.1;

    // 0D solver parameters. 
    // These are read in from the input JSON solver configuration file.
    double time_step_size_ = 0.0;
    int num_time_steps_ = 0;
    double absolute_tolerance_ = 0.0;
    int max_nliter_ = 0;
    int time_step_ = 0.0;
    int save_interval_counter_ = 0;
    int output_interval_ = 0;
    int system_size_ = 0;
    int num_output_steps_ = 0;
    int pts_per_cycle_ = 0;
    bool output_last_cycle_only_ = false;

    MODEL::Model<T>* model_;
    ALGEBRA::Integrator<T> integrator_;

    ALGEBRA::State<double> state_;
    std::vector<double> times_;
    std::vector<ALGEBRA::State<T>> states_;
};

