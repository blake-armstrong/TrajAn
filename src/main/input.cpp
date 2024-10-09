#include <CLI/CLI.hpp>
#include <memory>
#include <stdexcept>
#include <trajan/core/log.h>
#include <trajan/io/file.h>
#include <trajan/main/input.h>
#include <trajan/main/trajan.h>

namespace trajan::main {

namespace fs = std::filesystem;
// namespace core = trajan::core;
namespace io = trajan::io;

using FileHandler = std::unique_ptr<io::FileHandler>;
using TRAJAN = std::shared_ptr<trajan::TRAJAN>;

FileHandler read_input_file(const std::string &filename) {
  fs::path path = fs::path(filename);
  std::string ext = path.extension().string();
  trajan::log::debug("Attempting to read input from {}, file extension = {}",
                    filename, ext);
  if (!fs::exists(path)) {
    throw std::runtime_error(
        fmt::format("Input file does not exist: '{}'", filename));
  };

  FileHandler handler = io::get_handler(ext);
  handler->file_path = path;

  return handler;
}

void run_input_subcommand(InputOpts const &opt, TRAJAN &trajan) {
  std::vector<FileHandler> handlers;
  for (std::string file : opt.files) {
    trajan::log::info("File {} input.", file);
    handlers.push_back(read_input_file(file));
  }
  io::FileCompatibilityChecker checker;
  checker.check_and_modify_handlers(handlers);
}

CLI::App *add_input_subcommand(CLI::App &app, TRAJAN &trajan) {
  CLI::App *inp =
      app.add_subcommand("input", "Input coordinates and/or trajectory(s)");
  auto opts = std::make_shared<InputOpts>();
  inp->add_option("-f,--file", opts->files, "File name")->required();

  inp->callback([opts, &trajan]() { run_input_subcommand(*opts, trajan); });
  return inp;
}

} // namespace trajan::main
