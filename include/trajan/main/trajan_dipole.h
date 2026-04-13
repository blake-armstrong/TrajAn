#pragma once
#include <CLI/CLI.hpp>
#include <string>
#include <vector>
#include <trajan/core/pipeline.h>
#include <trajan/core/trajectory.h>

namespace trajan::main {

using trajan::core::Pipeline;
using trajan::core::Trajectory;

struct DipoleOpts {
  // Existing: whole-cell dipole output file
  std::string outfile;

  // Grid polarization field options
  std::string sel;              // molecule selector (required when grid_out set)
  std::string grid_out;         // output prefix for DX files
  std::vector<int> grid;        // explicit Na Nb Nc (overrides spacing)
  double spacing{1.0};          // target grid spacing in Å
  int order{4};                 // B-spline assignment order
};

void run_dipole_subcommand(const DipoleOpts &opts, Trajectory &traj,
                           const Pipeline &pipeline);
CLI::App *add_dipole_subcommand(CLI::App &app, Trajectory &traj,
                                Pipeline &pipeline);

} // namespace trajan::main
