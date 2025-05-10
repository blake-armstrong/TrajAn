#include <CLI/CLI.hpp>
#include <ankerl/unordered_dense.h>
#include <bitset>
#include <memory>
#include <trajan/core/log.h>
#include <trajan/core/neigh.h>
#include <trajan/core/util.h>
#include <trajan/io/file_handler.h>
#include <trajan/io/selection.h>
#include <trajan/io/text.h>
#include <trajan/main/trajan_rdf.h>

namespace trajan::main {

namespace core = trajan::core;
namespace io = trajan::io;
using uFileHandler = std::unique_ptr<io::FileHandler>;

void run_rdf_subcommand(const RDFOpts &opts) {
  std::vector<uFileHandler> handlers = io::read_input_files(opts.infiles);
  core::Trajectory trajectory;
  core::Frame &frame = trajectory.frame();

  std::vector<double> r(opts.nbins, 0.0);
  std::vector<size_t> nofr(opts.nbins, 0);
  std::vector<double> gofr(opts.nbins, 0.0);
  double bin_width = opts.rcut / opts.nbins;
  double inv_bin_width = 1 / bin_width;
  double norm = 4.0 * trajan::units::PI / 3.0;

  for (uFileHandler &handler : handlers) {
    bool read = handler->initialise();
    while (handler->read_frame(frame)) {
      core::UnitCell uc = frame.unit_cell();
      core::NeighbourList nl(uc, opts.rcut, opts.num_threads);
      io::SelectionCriteria sel1 = *opts.parsed_sel1;
      io::SelectionCriteria sel2 = *opts.parsed_sel2;
      std::vector<core::Entities> all_entities = {
          trajectory.get_entities(sel1), trajectory.get_entities(sel2)};
      trajan::log::debug("Entities found before deduplication = {}",
                         all_entities[0].size() + all_entities[1].size());
      auto result = trajan::util::combine_deduplicate_map(
          all_entities, core::VariantHash(), core::VariantEqual());
      core::Entities deduplicated_entities = result.first;
      std::vector<std::bitset<8>> presence_tracker = result.second;
      trajan::log::debug("Entities found after deduplication = {}",
                         deduplicated_entities.size());
      core::NeighbourListPacket nlp =
          trajectory.get_neighpack_from_entities(deduplicated_entities);
      nl.update(nlp);
      std::atomic<size_t> counter = 0;
      core::NeighbourCallback func =
          [&nofr, inv_bin_width, &presence_tracker, &counter](
              const core::Entity &ent1, const core::Entity &ent2, double rsq) {
            std::bitset<8> source1 = presence_tracker[ent1.idx];
            std::bitset<8> source2 = presence_tracker[ent2.idx];
            bool is_cross_section = (source1.test(0) && source2.test(1)) ||
                                    (source1.test(1) && source2.test(0));
            if (!is_cross_section) {
              return;
            }
            // don't need std::visit for typing
            // as it's inconsequential to rdf
            double r = std::sqrt(rsq);
            // trajan::log::debug("look {} {} {}", ent1.idx, ent2.idx, r);
            size_t bin_idx = r * inv_bin_width;
            nofr[bin_idx]++;
            counter++;
          };
      nl.iterate_neighbours(func);
      trajan::log::debug("Counter {}", static_cast<size_t>(counter));
      double sel2_density = all_entities[1].size() / uc.volume();
      double sel1_ref = all_entities[0].size() / 2.0;
      double density_norm = sel1_ref * sel2_density;
      for (size_t i = 0; i < opts.nbins; i++) {
        double ri = (i + 0.5) * bin_width;
        r[i] = ri;
        double shell_volume = norm * (std::pow(ri + bin_width / 2, 3) -
                                      std::pow(ri - bin_width / 2, 3));
        gofr[i] = nofr[i] / shell_volume / density_norm;
      }
    }
    handler->finalise();
  }
  TextFileWriter outfile;
  outfile.open(opts.outfile);
  outfile.write_line("{:>16} {:>16} {:>16}", "r", "nofr", "gofr");
  std::string fmt_str = "{:>16.8f} {:>16d} {:>16.8f}";
  for (size_t i = 0; i < opts.nbins; i++) {
    outfile.write_line(fmt_str, r[i], nofr[i], gofr[i]);
  }
  outfile.close();
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
