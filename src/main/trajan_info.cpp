#include <trajan/main/trajan_info.h>

namespace trajan::main {

void run_info_subcommand(const InfoOpts &opts, Trajectory &traj) {

  while (traj.next_frame()) {
    const auto uc = traj.unit_cell();
    if (uc) {
      size_t width = 24;
      const auto &ucv = uc.value();
      trajan::log::info("Unit cell information:");
      trajan::log::info("{:>{}}: {:>10.4f} {:>10.4f} {:>10.4f}",
                        "cell lengths (a b c)", width, ucv.a(), ucv.b(),
                        ucv.c());
      using occ::units::degrees;
      trajan::log::info("{:>{}}: {:>10.4f} {:>10.4f} {:>10.4f}",
                        "cell angles (α β γ)", width, degrees(ucv.alpha()),
                        degrees(ucv.beta()), degrees(ucv.gamma()));
      std::array<std::string, 3> labels{"A", "B", "C"};
      const auto &dt = ucv.direct();
      for (size_t i = 0; i < 3; i++) {
        const auto v = dt.col(i);
        trajan::log::info("{:>{}}: {:>10.4f} {:>10.4f} {:>10.4f}",
                          fmt::format("direct matrix ({})", labels[i]), width,
                          v(0), v(1), v(2));
      }
      const auto &rt = ucv.reciprocal();
      for (size_t i = 0; i < 3; i++) {
        const auto v = rt.col(i);
        const auto l = fmt::format("{}*", labels[i]);
        trajan::log::info("{:>{}}: {:>10.4f} {:>10.4f} {:>10.4f}",
                          fmt::format("reciprocal matrix ({})", l), width, v(0),
                          v(1), v(2));
      }
      trajan::log::info("{:>{}}: {}", "cell type", width, ucv.cell_type());
      trajan::log::info("{:>{}}: {:>14.4f}", "volume", width, ucv.volume());
      double area_bc = ucv.b() * ucv.c() * sin(ucv.alpha());
      double area_ac = ucv.a() * ucv.c() * sin(ucv.beta());
      double area_ab = ucv.a() * ucv.b() * sin(ucv.gamma());
      double surface_area = 2.0 * (area_bc + area_ac + area_ab);
      trajan::log::info("{:>{}}: {:>14.4f}", "surface area", width,
                        surface_area);
    } else {
      trajan::log::info("No unit cell information");
    }
  }
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
