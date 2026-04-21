#pragma once
#include <CLI/CLI.hpp>
#include <filesystem>
#include <optional>
#include <trajan/core/pipeline.h>
#include <trajan/core/trajectory.h>

namespace trajan::main {

using trajan::core::Pipeline;
using trajan::core::Trajectory;

struct SurfaceOpts {
  std::vector<int> hkl{0, 0, 1};
  std::optional<double> shift;
  double depth{1.0};
  bool atomic{false};
  std::vector<int> x_direction;
  std::filesystem::path outfile;
};

void run_surface_subcommand(const SurfaceOpts &opts, Trajectory &traj);
CLI::App *add_surface_subcommand(CLI::App &app, Trajectory &traj,
                                 Pipeline &pipeline);

} // namespace trajan::main
