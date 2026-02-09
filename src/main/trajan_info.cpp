#include <trajan/main/trajan_info.h>

namespace trajan::main {

void run_info_subcommand(const InfoOpts &opts, Trajectory &traj) {

};

CLI::App *add_info_subcommand(CLI::App &app, Trajectory &traj) {
  auto *info =
      app.add_subcommand("info", "Prints out a variety of information within "
                                 "the current Trajectory object.");
  auto opts = std::make_shared<InfoOpts>();
  info->add_option("--density", opts->density,
                   "Whether to output density info.");
  info->add_option("--cell", opts->cell, "Whether to output cell info.");
  info->add_option("--timings", opts->timings,
                   "Whether to output timing info.");

  info->callback([opts, &traj]() {
    trajan::log::set_subcommand_log_pattern("info");
    run_info_subcommand(*opts, traj);
  });

  return info;
};

} // namespace trajan::main
