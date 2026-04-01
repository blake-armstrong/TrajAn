#pragma once
#include <CLI/CLI.hpp>
#include <string>
#include <trajan/core/pipeline.h>
#include <trajan/core/trajectory.h>
#include <vector>

namespace trajan::main {

using trajan::core::Pipeline;
using trajan::core::Trajectory;

struct PmepotOpts {
  std::string outfile;
  double spacing{1.0};      // target grid spacing in Å (used if grid not set)
  std::vector<int> grid{};  // explicit Na Nb Nc (overrides spacing)
  double ewaldcof{0.25};    // Ewald coefficient β in Å⁻¹
  int order{4};             // B-spline order
};

void run_pmepot_subcommand(const PmepotOpts &opts, Trajectory &traj,
                           const Pipeline &pipeline);
CLI::App *add_pmepot_subcommand(CLI::App &app, Trajectory &traj,
                                Pipeline &pipeline);

} // namespace trajan::main
