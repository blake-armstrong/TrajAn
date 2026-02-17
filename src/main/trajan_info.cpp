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
    using ankerl::unordered_dense::map;
    map<std::string, size_t> atom_type_count;
    map<std::string, std::string> atom_type_elements;
    double mass = 0.0;
    size_t total_atoms = traj.atoms().size();
    for (const auto &a : traj.atoms()) {
      atom_type_count[a.type]++;
      atom_type_elements[a.type] = a.element.symbol();
      mass += a.element.mass();
    }
    trajan::log::info("Total atom count: {}", total_atoms);
    trajan::log::info("Unique atom types:");
    for (const auto &[at, count] : atom_type_count) {
      const auto &e = atom_type_elements[at];
      trajan::log::info("{:>5} ({:>2}): {:>10}", at, e, count);
    }
    size_t width = 12;
    trajan::log::info("Atom properties:");
    trajan::log::info("{:>{}}: {:>10.4f} g/mol", "total mass", width, mass);
    if (uc) {
      trajan::log::info("{:>{}}: {:>10.4f} g/cm^3", "density", width,
                        mass / uc->volume() * 1.6605388);
      trajan::log::info("{:>{}}: {:>10.4f} atom/Å^3 ", "density", width,
                        total_atoms / uc->volume());
    }
    core::Topology &top = traj.get_topology();
    top.print_summary();
    if (opts.detailed_top) {
      traj.frame().populate_angles(top);
      top.print_detailed();
    }
  }
};

CLI::App *add_info_subcommand(CLI::App &app, Trajectory &traj) {
  auto *info =
      app.add_subcommand("info", "Prints out a variety of information within "
                                 "the current Trajectory object.");
  auto opts = std::make_shared<InfoOpts>();
  info->add_flag("--detailed-topology", opts->detailed_top,
                 "Whether to output detailed topology info (if there is any "
                 "topology info).");
  info->add_option("--timings", opts->timings,
                   "Whether to output timing info.");

  info->callback([opts, &traj]() {
    trajan::log::set_subcommand_log_pattern("info");
    run_info_subcommand(*opts, traj);
  });

  return info;
};

} // namespace trajan::main
