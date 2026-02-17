#include <trajan/main/trajan_qc.h>

namespace trajan::main {

void run_qc_subcommand(QCOpts const &opts, Trajectory &traj) {

};

CLI::App *add_qc_subcommand(CLI::App &app, Trajectory &traj) {
  CLI::App *qc =
      app.add_subcommand("qc", "Quaternion clustering of carbonates");
  auto opts = std::make_shared<QCOpts>();
  qc->add_option("--carbonates", opts->carbonate_selection,
                 "First selection (prefix: i=atom indices, a=atom types, "
                 "j=molecule indices, m=molecule types)\n"
                 "Examples:\n"
                 "  j1,3-5      (molecule indices 1,3,4,5)\n"
                 "  mM1,M2      (molecule types M1,M2)")
      ->required();

  qc->callback([opts, &traj]() {
    trajan::log::set_subcommand_log_pattern("qc");
    run_qc_subcommand(*opts, traj);
  });
  return qc;
}

} // namespace trajan::main
