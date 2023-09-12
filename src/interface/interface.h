
#include "algebra/Integrator.h"
#include "algebra/SparseSystem.h"
#include "algebra/State.h"
#include "helpers/debug.h"
#include "helpers/helpers.h"
#include "io/io.h"
#include "io/csvwriter.h"
#include "model/Model.h"

#include <map>
#include <string>
#include <vector>

using S = algebra::SparseSystem;

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

    std::shared_ptr<zd_model::Model> model_;
    algebra::Integrator integrator_;

    algebra::State state_;
    std::vector<double> times_;
    std::vector<algebra::State> states_;
};

