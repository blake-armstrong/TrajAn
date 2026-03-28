#include <trajan/core/log.h>
#include <trajan/main/trajan_modify.h>

namespace trajan::main {

void register_modify_transforms(const ModifyOpts &opts, Pipeline &pipeline) {
  const double dx = opts.translate[0];
  const double dy = opts.translate[1];
  const double dz = opts.translate[2];

  if (dx != 0.0 || dy != 0.0 || dz != 0.0) {
    trajan::log::debug("modify: registering translate ({}, {}, {})", dx, dy,
                       dz);
    pipeline.add_transform([dx, dy, dz](trajan::core::Frame &frame) {
      for (auto &atom : frame.atoms()) {
        atom.x += dx;
        atom.y += dy;
        atom.z += dz;
      }
    });
  }

  if (opts.wrap) {
    trajan::log::debug("modify: registering PBC wrap");
    pipeline.add_transform([](trajan::core::Frame &frame) {
      if (!frame.has_unit_cell()) {
        trajan::log::warn(
            "modify --wrap: frame has no unit cell, skipping wrap");
        return;
      }
      const auto &uc = frame.unit_cell().value();
      for (auto &atom : frame.atoms()) {
        occ::Vec3 pos(atom.x, atom.y, atom.z);
        occ::Vec3 frac = uc.to_fractional(pos);
        frac = frac.array() - frac.array().floor();
        occ::Vec3 cart = uc.to_cartesian(frac);
        atom.x = cart[0];
        atom.y = cart[1];
        atom.z = cart[2];
      }
    });
  }
}

CLI::App *add_modify_subcommand(CLI::App &app, Trajectory &traj,
                                Pipeline &pipeline) {
  CLI::App *modify =
      app.add_subcommand("modify", "Register per-frame coordinate or topology "
                                   "modifications into the pipeline. Requires "
                                   "load subcommand.");

  auto opts = std::make_shared<ModifyOpts>();

  auto *trans_opt =
      modify->add_option("--translate", opts->translate,
                         "Translate all atom positions by dx dy dz (Angstroms)")
          ->expected(3);

  modify->add_flag("--wrap", opts->wrap,
                   "Wrap atom positions into the primary unit cell (requires "
                   "unit cell in trajectory)");

  modify->callback([opts, trans_opt, &pipeline]() {
    trajan::log::set_subcommand_log_pattern("modify");
    trajan::log::debug("Beginning modify subcommand");

    // If --translate was not given on the command line, keep defaults at zero.
    if (trans_opt->count() == 0) {
      opts->translate = {0.0, 0.0, 0.0};
    }

    register_modify_transforms(*opts, pipeline);
    trajan::log::debug("modify subcommand registered {} transform(s)",
                       pipeline.size());
  });

  // suppress "unused" warning — traj may be needed by future modify operations
  (void)traj;

  return modify;
}

} // namespace trajan::main
