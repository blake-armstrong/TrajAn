#include <CLI/CLI.hpp>

CLI::App *add_input_subcommand(CLI::App &app) {
  CLI::App *cg = app.add_subcommand("i", "Input coordinates and/or trajectory(s)");
  
}
