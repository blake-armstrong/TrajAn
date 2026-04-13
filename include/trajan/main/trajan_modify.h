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
  std::vector<std::string> raw_rotate{};             // [angle_deg, axis], empty if not set
  std::vector<double> rotate_origin{0.0, 0.0, 0.0}; // point to rotate about
  bool wrap{false};
  // Setter vectors — each entry pairs with the corresponding raw_sels entry.
  // If raw_sels is empty the single setter value applies to all atoms.
  std::vector<double> set_charges;
  std::vector<std::string> set_types;
  std::vector<std::string> set_mol_types;
  std::vector<std::string> raw_sels{};  // paired with set_charges/set_types/set_mol_types
  std::string raw_sel{};                // single sel for translate/rotate
  trajan::io::MolOrigin mol_origin{trajan::io::MolOrigin::CenterOfMass};
};

// Build and register all active transforms from opts into the pipeline.
// Called once at subcommand callback time (before the frame loop runs).
void register_modify_transforms(const ModifyOpts &opts, Trajectory &traj,
                                Pipeline &pipeline);

CLI::App *add_modify_subcommand(CLI::App &app, Trajectory &traj,
                                Pipeline &pipeline);

} // namespace trajan::main
