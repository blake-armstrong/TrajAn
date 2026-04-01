#pragma once
#include <CLI/CLI.hpp>
#include <optional>
#include <string>
#include <trajan/core/pipeline.h>
#include <trajan/core/trajectory.h>
#include <trajan/io/selection.h>
#include <vector>

namespace trajan::main {

using trajan::core::Pipeline;
using trajan::core::Trajectory;

struct ModifyOpts {
  std::vector<double> translate{0.0, 0.0, 0.0};
  std::vector<std::string> raw_rotate{};          // [angle_deg, axis], empty if not set
  std::vector<double> rotate_origin{0.0, 0.0, 0.0}; // point to rotate about
  bool wrap{false};
  std::optional<double> set_charge;
  std::optional<std::string> set_type;
  std::optional<std::string> set_mol_type;
  std::string raw_sel{};
  trajan::io::MolOrigin mol_origin{trajan::io::MolOrigin::CenterOfMass};
};

// Build and register all active transforms from opts into the pipeline.
// Called once at subcommand callback time (before the frame loop runs).
void register_modify_transforms(const ModifyOpts &opts, Trajectory &traj,
                                Pipeline &pipeline);

CLI::App *add_modify_subcommand(CLI::App &app, Trajectory &traj,
                                Pipeline &pipeline);

} // namespace trajan::main
