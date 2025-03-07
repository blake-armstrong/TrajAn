#include <CLI/CLI.hpp>
#include <memory>
#include <trajan/core/log.h>
#include <trajan/core/neigh.h>
#include <trajan/io/file.h>
#include <trajan/io/selection.h>
#include <trajan/main/trajan_rdf.h>

namespace trajan::main {

namespace core = trajan::core;
namespace io = trajan::io;
using uFileHandler = std::unique_ptr<io::FileHandler>;

void run_rdf_subcommand(RDFOpts const &opts) {
  std::vector<uFileHandler> handlers = io::read_input_files(opts.infiles);
  core::Frame frame;
  for (uFileHandler &handler : handlers) {
    while (handler->read_frame(frame)) {
      core::UnitCell uc = frame.get_uc();
      // TODO: make a way to estimate num_threads
      // depending on this system size and num atoms
      int num_threads = 1;
      core::NeighbourList neigh_list(uc, opts.rcut, num_threads);
      std::vector<core::Atom> atoms = frame.get_atoms();
      Mat3N atom_pos = frame.get_atom_positions();
      neigh_list.update(atoms, atom_pos);
      size_t Nb = 0;
      // for (core::Atom &atom : atoms) {
      //   if (!atom)
      // }
      // get unit cell volume

      // process_frame(frame);
    }
  }
}

CLI::App *add_rdf_subcommand(CLI::App &app) {
  CLI::App *rdf =
      app.add_subcommand("rdf", "Radial Pair Distribution Function");
  auto opts = std::make_shared<RDFOpts>();
  rdf->add_option("--tr,--traj", opts->infiles, "Input trajectory file name")
      ->required()
      ->check(CLI::ExistingFile);
  rdf->add_option("--o,--out", opts->outfile, "Output file for RDF data")
      ->capture_default_str();
  rdf->add_option("--rc,--rcut", opts->rcut, "RDF cutoff")
      ->capture_default_str();
  rdf->add_option("--nb,--nbins", opts->nbins, "Number of bins for RDF")
      ->capture_default_str();
  std::string sel1 = "--s1,--sel1";
  rdf->add_option(sel1, opts->raw_sel1,
                  "First selection (prefix: i=indices, t=types, m=molecules)\n"
                  "Examples:\n"
                  "  i1,2,3-5    (indices 1,2,3,4,5)\n"
                  "  tC,N,O      (atom types C, N, O)\n"
                  "  m1,3-5      (molecules 1,3,4,5)")
      ->required()
      ->check(io::selection_validator(opts->parsed_sel1));

  rdf->add_option("--s2,--sel2", opts->raw_sel2,
                  fmt::format("Second selection (same format as {})", sel1))
      ->required()
      ->check(io::selection_validator(opts->parsed_sel2));

  rdf->callback([opts]() { run_rdf_subcommand(*opts); });
  return rdf;
}

} // namespace trajan::main
