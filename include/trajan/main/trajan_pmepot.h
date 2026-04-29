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
  bool real_space{false};   // use direct real-space Coulomb sum instead of PME
  double cutoff{10.0};      // neighbour cutoff for real-space sum (Å)
  bool grid_spread{false};  // pmepot_custom style: B-spline charge grid + grid-space sum
  int cell_radius{2};       // ±cell search radius for --grid-spread
};

void run_pmepot_subcommand(const PmepotOpts &opts, Trajectory &traj,
                           const Pipeline &pipeline);
CLI::App *add_pmepot_subcommand(CLI::App &app, Trajectory &traj,
                                Pipeline &pipeline);

struct DXReduceOpts {
  std::string infile;
  std::string outfile;
  std::string average_axes; // e.g. "a", "bc", "ab" — axes to collapse by averaging
};

void run_dxreduce_subcommand(const DXReduceOpts &opts);
CLI::App *add_dxreduce_subcommand(CLI::App &app);

} // namespace trajan::main
