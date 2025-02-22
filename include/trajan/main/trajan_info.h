#pragma once

#include <CLI/CLI.hpp>
#include <string>
#include <vector>

namespace trajan::main {
namespace fs = std::filesystem;
struct InfoOpts {
  fs::path infile;
};

void run_info_subcommand(const InfoOpts &opts);
CLI::App *add_info_subcommand(CLI::App &app);

} // namespace trajan::main
