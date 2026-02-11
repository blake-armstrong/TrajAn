#pragma once
#include <CLI/CLI.hpp>
#include <trajan/core/trajectory.h>

namespace trajan::main {

using trajan::core::Trajectory;

struct QCOpts {
  std::string carbonate_selection;
};

void run_qc_subcommand(QCOpts const &opts, Trajectory &traj);
CLI::App *add_qc_subcommand(CLI::App &app, Trajectory &traj);

} // namespace trajan::main
