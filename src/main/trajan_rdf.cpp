#include <CLI/CLI.hpp>
#include <memory>
#include <trajan/core/log.h>
#include <trajan/core/neigh.h>
#include <trajan/core/util.h>
#include <trajan/io/file.h>
#include <trajan/io/selection.h>
#include <trajan/main/trajan_rdf.h>
#include <variant>

namespace trajan::main {

namespace core = trajan::core;
namespace io = trajan::io;
using uFileHandler = std::unique_ptr<io::FileHandler>;

void run_rdf_subcommand(const RDFOpts &opts) {
  std::vector<uFileHandler> handlers = io::read_input_files(opts.infiles);
  core::Trajectory trajectory;
  core::Frame frame = trajectory.frame();
  for (uFileHandler &handler : handlers) {
    bool read = handler->initialise();
    if (!read) {
      throw std::runtime_error(fmt::format(
          "Unable to open file for reading: '{}'", handler->file_name()));
    }
    trajan::log::debug(
        fmt::format("Successfully opened file '{}'", handler->file_name()));

    while (handler->read_frame(frame)) {
      trajan::log::debug(fmt::format("uc: {}", frame.unit_cell().dummy()));
      core::NeighbourList nl(frame.unit_cell(), opts.rcut, opts.num_threads);
      io::SelectionCriteria sel1 = *opts.parsed_sel1;
      core::Entities ents1 = trajectory.get_entities(sel1);
      io::SelectionCriteria sel2 = *opts.parsed_sel2;
      core::Entities ents2 = trajectory.get_entities(sel2);
      size_t idx_cutoff = ents1.size();
      core::Entities combined_ents =
          trajan::util::combine_vectors<core::EntityType>({ents1, ents2});
      core::NeighbourListPacket nlp =
          trajectory.get_neighpack_from_entities(combined_ents);
      nl.update(nlp);
      std::vector<double> rdf_hist(opts.nbins, 0.0);
      core::NeighbourCallback func = [&combined_ents, idx_cutoff, &rdf_hist](
                                         const core::Entity &ent1,
                                         const core::Entity &ent2, double rsq) {
        bool is_cross_selection =
            (ent1.idx < idx_cutoff && ent2.idx >= idx_cutoff) ||
            (ent1.idx >= idx_cutoff && ent2.idx < idx_cutoff);
        if (!is_cross_selection) {
          return;
        }
        // const auto &sel1_ent = (ent1.idx < idx_cutoff)
        //                            ? combined_ents[ent1.idx]
        //                            : combined_ents[ent2.idx];
        //
        // const auto &sel2_ent = (ent1.idx >= idx_cutoff)
        //                            ? combined_ents[ent1.idx]
        //                            : combined_ents[ent2.idx];
        // don't need std::visit for typing
        // as it's inconsequential to rdf
        double r = std::sqrt(rsq);
        // TODO: do something lmao
      };
    }
    handler->finalise();
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

  opts->num_threads = app.get_option("--threads")->as<size_t>();

  rdf->callback([opts]() { run_rdf_subcommand(*opts); });
  return rdf;
}

} // namespace trajan::main
