#pragma once
#include <CLI/CLI.hpp>
#include <filesystem>
#include <trajan/core/pipeline.h>
#include <trajan/core/trajectory.h>

namespace trajan::main {

using trajan::core::Pipeline;
using trajan::core::Trajectory;

struct WriteOpts {
  std::filesystem::path outfile;
  bool original_ids{false};
};

void run_write_subcommand(const WriteOpts &opts, Trajectory &traj,
                          const Pipeline &pipeline);
CLI::App *add_write_subcommand(CLI::App &app, Trajectory &traj,
                               Pipeline &pipeline);

} // namespace trajan::main
