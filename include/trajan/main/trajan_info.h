#pragma once

#include <CLI/CLI.hpp>
#include <trajan/core/trajectory.h>

namespace trajan::main {

using trajan::core::Trajectory;

struct InfoOpts {
  bool density{true};
  bool cell{true};
  bool timings{true};
};

void run_info_subcommand(const InfoOpts &opts, Trajectory &traj);
CLI::App *add_info_subcommand(CLI::App &app, Trajectory &traj);

} // namespace trajan::main
