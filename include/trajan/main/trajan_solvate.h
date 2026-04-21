#pragma once
#include <CLI/CLI.hpp>
#include <filesystem>
#include <optional>
#include <trajan/core/pipeline.h>
#include <trajan/core/trajectory.h>

namespace trajan::main {

using trajan::core::Pipeline;
using trajan::core::Trajectory;

struct SolvateOpts {
  std::filesystem::path molecule_file;
  std::vector<double> origin{0.0, 0.0, 0.0};
  std::vector<double> box{10.0, 10.0, 10.0};
  double min_dist{2.5};
  std::optional<int> n_molecules;
  std::optional<double> density;
  int max_attempts{200};
  unsigned int seed{42};
};

void run_solvate_subcommand(const SolvateOpts &opts, Trajectory &traj,
                            Pipeline &pipeline);
CLI::App *add_solvate_subcommand(CLI::App &app, Trajectory &traj,
                                 Pipeline &pipeline);

} // namespace trajan::main
