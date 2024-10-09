#pragma once

#include <CLI/CLI.hpp>
#include <memory>
#include <string>
#include <trajan/core/trajectory.h>
#include <trajan/io/file.h>
#include <trajan/main/trajan.h>
#include <vector>

namespace trajan::main {

namespace core = trajan::core;
namespace io = trajan::io;

struct InputOpts {
  std::vector<std::string> files;
};

std::unique_ptr<io::FileHandler> read_input_file(const std::string &filename);
// trajan::core::Trajectory
// combine_trajectories(std::vector<trajan::core::Trajectory> trajectories);
void run_input_subcommand(InputOpts const &opt,
                          std::shared_ptr<trajan::TRAJAN> &trajan);
CLI::App *add_input_subcommand(CLI::App &app,
                               std::shared_ptr<trajan::TRAJAN> &trajan);

} // namespace trajan::main
