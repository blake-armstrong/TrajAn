#include <CLI/CLI.hpp>
#include <tacpp/core/log.h>
#include <tacpp/main/tacpp_input.h>

namespace tacpp::main {

CLI::App *add_input_subcommand(CLI::App &app) {
  CLI::App *inp =
      app.add_subcommand("input", "Input coordinates and/or trajectory(s)");
  auto opts = std::make_shared<InputOpts>();
  inp->add_option("-f,--file", opts->files, "File name")->required();

  inp->callback([opts]() { run_input_subcommand(*opts); });
  return inp;
}

void run_input_subcommand(InputOpts const &opt) {
  for (std::string file : opt.files) {
    tacpp::log::info("File {} input.", file);
  }
}

} // namespace tacpp::main
