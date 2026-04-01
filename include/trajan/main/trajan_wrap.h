#pragma once
#include <CLI/CLI.hpp>
#include <trajan/core/pipeline.h>
#include <trajan/core/trajectory.h>

namespace trajan::main {

using trajan::core::Pipeline;
using trajan::core::Trajectory;

struct WrapOpts {
  bool atomic{false};
};

void register_wrap_transform(const WrapOpts &opts, Trajectory &traj,
                             Pipeline &pipeline);
CLI::App *add_wrap_subcommand(CLI::App &app, Trajectory &traj,
                              Pipeline &pipeline);

} // namespace trajan::main
