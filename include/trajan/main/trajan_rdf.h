#pragma once

#include <CLI/CLI.hpp>
#include <string>
#include <trajan/io/file.h>
#include <trajan/io/selection.h>
#include <vector>

namespace trajan::main {

// namespace core = trajan::core;
namespace io = trajan::io;
namespace fs = std::filesystem;

struct RDFOpts {
  std::vector<fs::path> infiles;
  std::string outfile = "gofr.out";
  double rcut = 6.0;
  int nbins = 100;
  std::string raw_sel1, raw_sel2;
  std::optional<io::SelectionCriteria> parsed_sel1, parsed_sel2;
};

void run_rdf_subcommand(RDFOpts const &opts);
CLI::App *add_rdf_subcommand(CLI::App &app);

} // namespace trajan::main
