#include <CLI/CLI.hpp>
#include <trajan/io/file.h>
#include <trajan/main/trajan_info.h>

namespace trajan::main {
namespace io = trajan::io;
using uFileHandler = std::unique_ptr<io::FileHandler>;

void run_info_subcommand(const InfoOpts &opts) {
  uFileHandler handler = io::read_input_file(opts.infile);
  core::Frame frame;
  handler->read_frame(frame);
}

CLI::App *add_info_subcommand(CLI::App &app) {
  CLI::App *info = app.add_subcommand("info", "General System Information");
  auto opts = std::make_shared<InfoOpts>();
  info->add_option("--f,--file", opts->infile, "Input file name")
      ->required()
      ->check(CLI::ExistingFile);
  info->callback([opts]() { run_info_subcommand(*opts); });
  return info;
}

} // namespace trajan::main
