// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "Model.h"

template <typename block_type>
BlockFactoryFunc block_factory() {
  return [](int count, Model* model) -> Block* {
    return new block_type(count, model);
  };
}

Model::Model() {
  // Add all implemented blocks to factory
  block_factory_map = {
      {"BloodVessel", block_factory<BloodVessel>()},
      {"ChamberSphere", block_factory<ChamberSphere>()},
      {"BloodVesselJunction", block_factory<BloodVesselJunction>()},
      {"ClosedLoopCoronaryLeft", block_factory<ClosedLoopCoronaryLeftBC>()},
      {"ClosedLoopCoronaryRight", block_factory<ClosedLoopCoronaryRightBC>()},
      {"ClosedLoopHeartAndPulmonary",
       block_factory<ClosedLoopHeartPulmonary>()},
      {"ClosedLoopRCR", block_factory<ClosedLoopRCRBC>()},
      {"CORONARY", block_factory<OpenLoopCoronaryBC>()},
      {"CORONARY_VAR_RES", block_factory<OpenLoopCoronaryVarResBC>()},
      {"FLOW", block_factory<FlowReferenceBC>()},
      {"IMPEDANCE", block_factory<ImpedanceBC>()},
      {"NORMAL_JUNCTION", block_factory<Junction>()},
      {"PRESSURE", block_factory<PressureReferenceBC>()},
      {"RCR", block_factory<WindkesselBC>()},
      {"RESISTANCE", block_factory<ResistanceBC>()},
      {"resistive_junction", block_factory<ResistiveJunction>()},
      {"ValveTanh", block_factory<ValveTanh>()},
      {"ChamberElastanceInductor", block_factory<ChamberElastanceInductor>()},
      {"ChamberElastanceInductorExponential",
       block_factory<ChamberElastanceInductorExponential>()},
      {"BloodVesselCRL", block_factory<BloodVesselCRL>()},
      {"PiecewiseValve", block_factory<PiecewiseValve>()},
      {"LinearElastanceChamber", block_factory<LinearElastanceChamber>()},
      {"BloodVesselRC", block_factory<BloodVesselRC>()}};
}

Model::~Model() {}

Block* Model::create_block(const std::string& block_type) {
  // Get block from factory
  auto it = block_factory_map.find(block_type);
  if (it == block_factory_map.end()) {
    throw std::runtime_error("Invalid block type " + block_type);
  }
  Block* block = it->second(block_count, this);
  return block;
}

int Model::add_block(Block* block, const std::string_view& name,
                     const std::vector<int>& block_param_ids, bool internal) {
  // Set global parameter IDs
  block->setup_params_(block_param_ids);

  auto name_string = static_cast<std::string>(name);

  if (internal) {
    hidden_blocks.push_back(std::shared_ptr<Block>(block));
  } else {
    blocks.push_back(std::shared_ptr<Block>(block));
  }

  block_types.push_back(block->block_type);
  block_index_map.insert({name_string, block_count});
  block_names.push_back(name_string);

  return block_count++;
}

int Model::add_block(const std::string& block_name,
                     const std::vector<int>& block_param_ids,
                     const std::string_view& name, bool internal) {
  // Generate block from factory
  auto block = this->create_block(block_name);

  // Add block to model
  return this->add_block(block, name, block_param_ids, internal);
}

bool Model::has_block(const std::string& name) const {
  if (block_index_map.find(name) == block_index_map.end()) {
    return false;
  } else {
    return true;
  }
}

Block* Model::get_block(const std::string_view& name) const {
  auto name_string = static_cast<std::string>(name);

  if (!has_block(name_string)) {
    throw std::runtime_error("No block defined with name " + name_string);
  }

  return blocks[block_index_map.at(name_string)].get();
}

Block* Model::get_block(int block_id) const {
  if (block_id >= blocks.size()) {
    return hidden_blocks[block_id - blocks.size()].get();
  }

  return blocks[block_id].get();
}

BlockType Model::get_block_type(const std::string_view& name) const {
  auto name_string = static_cast<std::string>(name);

  if (block_index_map.find(name_string) == block_index_map.end()) {
    throw std::runtime_error("Could not find block with name " + name_string);
  }

  return block_types[block_index_map.at(name_string)];
}

std::string Model::get_block_name(int block_id) const {
  return block_names[block_id];
}

int Model::add_node(const std::vector<Block*>& inlet_eles,
                    const std::vector<Block*>& outlet_eles,
                    const std::string_view& name) {
  DEBUG_MSG("Adding node " << name);
  auto node = std::shared_ptr<Node>(
      new Node(node_count, inlet_eles, outlet_eles, this));
  nodes.push_back(node);
  node_names.push_back(static_cast<std::string>(name));

  return node_count++;
}

std::string Model::get_node_name(int node_id) const {
  return node_names[node_id];
}

int Model::add_parameter(double value) {
  parameters.push_back(Parameter(parameter_count, value));
  parameter_values.push_back(parameters.back().get(0.0));
  return parameter_count++;
}

int Model::add_parameter(const std::vector<double>& times,
                         const std::vector<double>& values, bool periodic) {
  auto param = Parameter(parameter_count, times, values, periodic);
  if (periodic && (param.is_constant == false)) {
    if ((this->cardiac_cycle_period > 0.0) &&
        (param.cycle_period != this->cardiac_cycle_period)) {
      throw std::runtime_error(
          "Inconsistent cardiac cycle period defined in parameters");
    }
    this->cardiac_cycle_period = param.cycle_period;
  }
  parameter_values.push_back(param.get(0.0));
  parameters.push_back(std::move(param));
  return parameter_count++;
}

Parameter* Model::get_parameter(int param_id) { return &parameters[param_id]; }

double Model::get_parameter_value(int param_id) const {
  return parameter_values[param_id];
}

void Model::update_parameter_value(int param_id, double param_value) {
  parameter_values[param_id] = param_value;
}

namespace {

/// Upper connection limit meaning "any number".
constexpr int kUnlimited = -1;

/**
 * @brief Inlet and outlet connections a block type needs
 *
 * Blocks address the pressure and flow of their nodes by position in
 * `global_var_ids` (see Block::setup_dofs_), so these are the connection
 * counts that each block's setup_dofs() and update_constant() assume.
 */
struct ConnectivityRule {
  int min_inlets;   ///< Smallest valid number of inlet connections
  int max_inlets;   ///< Largest valid number, or kUnlimited
  int min_outlets;  ///< Smallest valid number of outlet connections
  int max_outlets;  ///< Largest valid number, or kUnlimited
  /// Block terminates the model: exactly one connection, on either side
  bool single_sided{false};
};

/**
 * @brief Connections required by a block type
 *
 * Written as a switch over all block types without a `default` case, so a
 * newly added block type shows up as a compiler warning here instead of as a
 * missing connectivity check.
 *
 * @param type Type of the block
 * @return Connections the block type requires
 */
ConnectivityRule get_connectivity_rule(BlockType type) {
  switch (type) {
    // Vessels, valves and chambers sit between exactly two nodes.
    case BlockType::blood_vessel:
    case BlockType::blood_vessel_CRL:
    case BlockType::blood_vessel_rc:
    case BlockType::valve_tanh:
    case BlockType::piecewise_valve:
    case BlockType::chamber_elastance_inductor:
    case BlockType::chamber_elastance_inductor_exponential:
    case BlockType::chamber_sphere:
    case BlockType::linear_elastance_chamber:
    // Closed-loop blocks sit inside the circulation, so they have both an
    // upstream and a downstream connection.
    case BlockType::closed_loop_rcr_bc:
    case BlockType::closed_loop_coronary_left_bc:
    case BlockType::closed_loop_coronary_right_bc:
    case BlockType::closed_loop_heart_pulmonary:
      return {1, 1, 1, 1};

    // Junctions distribute flow between any number of inlets and outlets.
    case BlockType::junction:
    case BlockType::resistive_junction:
      return {1, kUnlimited, 1, kUnlimited};

    // A blood vessel junction models each of its outlets as a vessel of its
    // own, all fed by a single inlet.
    case BlockType::blood_vessel_junction:
      return {1, 1, 1, kUnlimited};

    // Open-loop boundary conditions terminate the model on one side.
    case BlockType::flow_bc:
    case BlockType::pressure_bc:
    case BlockType::resistance_bc:
    case BlockType::windkessel_bc:
    case BlockType::impedance_bc:
    case BlockType::open_loop_coronary_bc:
    case BlockType::open_loop_coronary_var_res_bc:
      return {0, 1, 0, 1, true};
  }
  throw std::runtime_error(
      "Internal error: no connectivity rule for block type " +
      std::to_string(static_cast<int>(type)) + ".");
}

/**
 * @brief Name under which a block type is registered in the block factory
 *
 * Read from the factory so that error messages name block types the way the
 * configuration file does. One throwaway block is constructed per registered
 * type; block constructors only initialize members, so this doesn't touch
 * `model`.
 *
 * @param model Model whose block factory is queried
 * @param type Type of the block
 * @return Registered type name, or a placeholder if the type is unregistered
 */
std::string get_block_type_name(Model& model, BlockType type) {
  for (const auto& [type_name, factory] : model.block_factory_map) {
    std::unique_ptr<Block> block(factory(0, &model));
    if (block->block_type == type) {
      return std::string(type_name);
    }
  }
  return "unknown type";
}

/// Phrase a connection count limit, e.g. "exactly 1" or "at least 1".
std::string describe_limit(int min_count, int max_count) {
  if (min_count == max_count) {
    return "exactly " + std::to_string(min_count);
  }
  if (max_count == kUnlimited) {
    return "at least " + std::to_string(min_count);
  }
  return "between " + std::to_string(min_count) + " and " +
         std::to_string(max_count);
}

/// Phrase a number of connections, e.g. "1 inlet" or "2 inlets".
std::string describe_count(size_t count, const std::string& noun) {
  return std::to_string(count) + " " + noun + (count == 1 ? "" : "s");
}

}  // namespace

void Model::validate_connectivity() {
  std::vector<std::string> problems;

  // Check that every block has the connections its equations assume.
  for (const auto& block : blocks) {
    const ConnectivityRule rule = get_connectivity_rule(block->block_type);
    const size_t num_inlets = block->inlet_nodes.size();
    const size_t num_outlets = block->outlet_nodes.size();

    if (rule.single_sided) {
      // A boundary condition attaches with a single connection, on the side
      // of the block it terminates.
      if (num_inlets + num_outlets != 1) {
        problems.push_back(
            "Block '" + block->get_name() + "' (" +
            get_block_type_name(*this, block->block_type) + ") has " +
            describe_count(num_inlets + num_outlets, "connection") +
            ", but a boundary condition needs exactly 1, as the inlet or the "
            "outlet of one other block.");
      }
      continue;
    }

    const bool inlets_ok = static_cast<int>(num_inlets) >= rule.min_inlets &&
                           (rule.max_inlets == kUnlimited ||
                            static_cast<int>(num_inlets) <= rule.max_inlets);
    const bool outlets_ok = static_cast<int>(num_outlets) >= rule.min_outlets &&
                            (rule.max_outlets == kUnlimited ||
                             static_cast<int>(num_outlets) <= rule.max_outlets);
    if (inlets_ok && outlets_ok) {
      continue;
    }
    problems.push_back(
        "Block '" + block->get_name() + "' (" +
        get_block_type_name(*this, block->block_type) + ") has " +
        describe_count(num_inlets, "inlet connection") + " and " +
        describe_count(num_outlets, "outlet connection") + ", but needs " +
        describe_limit(rule.min_inlets, rule.max_inlets) + " inlet and " +
        describe_limit(rule.min_outlets, rule.max_outlets) + " outlet.");
  }

  if (problems.empty()) {
    return;
  }
  std::string message = "Invalid connectivity in the model:";
  for (const auto& problem : problems) {
    message += "\n  - " + problem;
  }
  throw std::runtime_error(message);
}

void Model::finalize() {
  DEBUG_MSG("Check connectivity of blocks");
  validate_connectivity();
  DEBUG_MSG("Setup degrees-of-freedom of nodes");
  for (auto& node : nodes) {
    node->setup_dofs(dofhandler);
  }
  DEBUG_MSG("Setup degrees-of-freedom of blocks");
  for (auto& block : blocks) {
    block->setup_dofs(dofhandler);
  }
  DEBUG_MSG("Setup model-dependent parameters");
  for (auto& block : blocks) {
    block->setup_model_dependent_params();
  }
}

int Model::get_num_blocks(bool internal) const {
  int num_blocks = blocks.size();

  if (internal) {
    num_blocks += hidden_blocks.size();
  }

  return num_blocks;
}

void Model::update_constant(SparseSystem& system) {
  for (auto block : blocks) {
    block->update_constant(system, parameter_values);
  }
}

void Model::update_time(SparseSystem& system, double time) {
  this->time = time;

  for (auto& param : parameters) {
    parameter_values[param.id] = param.get(time);
  }

  for (auto block : blocks) {
    block->update_time(system, parameter_values);
  }
}

void Model::update_solution(SparseSystem& system,
                            Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
                            Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) {
  for (auto block : blocks) {
    block->update_solution(system, parameter_values, y, dy);
  }
}

void Model::post_solve(Eigen::Matrix<double, Eigen::Dynamic, 1>& y) {
  for (auto block : blocks) {
    block->post_solve(y);
  }
}

void Model::accept_timestep(const Eigen::Matrix<double, Eigen::Dynamic, 1>& y) {
  for (size_t i = 0; i < get_num_blocks(true); i++) {
    get_block(i)->accept_timestep(y);
  }
}

std::vector<std::vector<double>> Model::get_persistent_states() const {
  std::vector<std::vector<double>> states;
  states.reserve(get_num_blocks(true));

  for (size_t i = 0; i < get_num_blocks(true); i++) {
    states.push_back(get_block(i)->get_persistent_state());
  }

  return states;
}

void Model::set_persistent_states(
    const std::vector<std::vector<double>>& states) {
  if (states.size() != static_cast<size_t>(get_num_blocks(true))) {
    throw std::runtime_error(
        "Persistent-state block count mismatch in "
        "Model::set_persistent_states");
  }

  for (size_t i = 0; i < states.size(); i++) {
    get_block(i)->set_persistent_state(states[i]);
  }
}

void Model::clear_persistent_states() {
  for (size_t i = 0; i < get_num_blocks(true); i++) {
    get_block(i)->set_persistent_state({});
  }
}

void Model::set_time_step_size(double time_step_size) {
  this->time_step_size = time_step_size;
}

double Model::get_time_step_size() const { return time_step_size; }

void Model::to_steady() {
  for (auto& param : parameters) {
    param.to_steady();
  }

  // Special handling for time-varying capacitance
  for (size_t i = 0; i < get_num_blocks(true); i++) {
    get_block(i)->steady = true;
    if ((block_types[i] == BlockType::windkessel_bc) ||
        (block_types[i] == BlockType::closed_loop_rcr_bc)) {
      int param_id_capacitance = blocks[i]->global_param_ids[1];
      double value = parameters[param_id_capacitance].get(0.0);
      param_value_cache.insert({param_id_capacitance, value});
      parameters[param_id_capacitance].update(0.0);
    }
  }
}

void Model::to_unsteady() {
  for (auto& param : parameters) {
    param.to_unsteady();
  }
  for (auto& [param_id_capacitance, value] : param_value_cache) {
    // DEBUG_MSG("Setting Windkessel capacitance back to " << value);
    parameters[param_id_capacitance].update(value);
  }
  for (size_t i = 0; i < get_num_blocks(true); i++) {
    get_block(i)->steady = false;
  }
}

TripletsContributions Model::get_num_triplets() const {
  TripletsContributions triplets_sum;

  for (auto& elem : blocks) {
    triplets_sum += elem->get_num_triplets();
  }

  return triplets_sum;
}

void Model::setup_initial_state_dependent_parameters(State initial_state) {
  DEBUG_MSG("Setup initial state dependent parameters");
  for (auto& block : blocks) {
    block->setup_initial_state_dependent_params(initial_state,
                                                parameter_values);
  }
}

void Model::update_has_windkessel_bc(bool has_windkessel) {
  has_windkessel_bc = has_windkessel;
}

void Model::update_largest_windkessel_time_constant(double time_constant) {
  largest_windkessel_time_constant = time_constant;
}

bool Model::get_has_windkessel_bc() { return has_windkessel_bc; }

double Model::get_largest_windkessel_time_constant() {
  return largest_windkessel_time_constant;
}
