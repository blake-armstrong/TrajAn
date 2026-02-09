#pragma once
#include <CLI/CLI.hpp>
#include <trajan/core/trajectory.h>
#include <trajan/core/util.h>
#include <trajan/io/selection.h>

namespace trajan::main {

using trajan::core::Trajectory;

struct BondCriteria {
  std::vector<io::SelectionCriteria> sel1;
  std::vector<io::SelectionCriteria> sel2;
  enum class ComparisonOp { LessThan, GreaterThan } op;
  double threshold;
};

struct TopologyOpts {
  bool top_auto{true};
  int update_frequency{0};
  std::vector<std::string> nb_raw_sel;
  std::vector<std::string> bc_raw_sel;
  // std::vector<io::SelectionCriteria> nb_parsed_sel;
  // std::vector<BondCriteria> bond_criterias;
  double bond_tolerance;
};

const auto MOLECULE_RESTRICTIONS =
    std::make_optional<std::vector<char>>({'j', 'm'});

BondCriteria bond_criteria_validator(
    const std::string &input,
    std::optional<std::vector<char>> restrictions = std::nullopt);

void run_topology_subcommand(const TopologyOpts &opts, Trajectory &traj);
CLI::App *add_topology_subcommand(CLI::App &app, Trajectory &traj);

} // namespace trajan::main
