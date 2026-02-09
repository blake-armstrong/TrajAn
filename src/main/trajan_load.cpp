#include <trajan/main/trajan_load.h>

namespace trajan::main {

void run_load_subcommand(const LoadOpts &opts, Trajectory &traj) {
  if (opts.into_mem) {
    trajan::log::debug("Files will be loaded into memory");
    traj.load_files_into_memory(opts.infiles);
    return;
  }
  trajan::log::debug("Files will NOT be loaded into memory.");
  trajan::log::debug("File paths saved for analysis in next subcommand.");
  traj.load_files(opts.infiles);
}

CLI::App *add_load_subcommand(CLI::App &app, Trajectory &traj) {
  CLI::App *load =
      app.add_subcommand("load", "Load trajectory data into program for "
                                 "analysis. Required by most subcommands.");
  auto opts = std::make_shared<LoadOpts>();
  load->add_option("files", opts->infiles, "Input trajectory file names")
      // ->required()
      ->check(CLI::ExistingFile);
  load->add_flag("--into-mem", opts->into_mem,
                 "Whether or not to load files into memory. This is not "
                 "recommended for command-line use.");
  load->callback([opts, &traj]() {
    trajan::log::set_subcommand_log_pattern("load");
    trajan::log::debug("Beginning load subcommand");
    run_load_subcommand(*opts, traj);
    trajan::log::debug("load subcommand completed successfully");
  });
  return load;
}

} // namespace trajan::main
