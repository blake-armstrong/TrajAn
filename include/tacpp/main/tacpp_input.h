#pragma once

#include <CLI/CLI.hpp>
#include <string>
#include <vector>

namespace tacpp::main {

struct InputOpts {
  std::vector<std::string> files;
};

CLI::App *add_input_subcommand(CLI::App &app);
void run_input_subcommand(InputOpts const &opt);

} // namespace tacpp::main
