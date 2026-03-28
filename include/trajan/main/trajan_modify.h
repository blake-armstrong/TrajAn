#pragma once
#include <CLI/CLI.hpp>
#include <trajan/core/pipeline.h>
#include <trajan/core/trajectory.h>
#include <vector>

namespace trajan::main {

using trajan::core::Pipeline;
using trajan::core::Trajectory;

struct ModifyOpts {
  std::vector<double> translate{0.0, 0.0, 0.0};
  bool wrap{false};
};

// Build and register all active transforms from opts into the pipeline.
// Called once at subcommand callback time (before the frame loop runs).
void register_modify_transforms(const ModifyOpts &opts, Pipeline &pipeline);

CLI::App *add_modify_subcommand(CLI::App &app, Trajectory &traj,
                                Pipeline &pipeline);

} // namespace trajan::main
