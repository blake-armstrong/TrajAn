#include <trajan/core/log.h>
#include <trajan/main/trajan_wrap.h>

namespace trajan::main {

void register_wrap_transform(const WrapOpts &opts, Trajectory &traj,
                             Pipeline &pipeline) {
  if (opts.atomic) {
    trajan::log::debug("wrap: registering atomic PBC wrap");
    pipeline.add_transform("wrap", [](trajan::core::Frame &frame) {
      if (!frame.has_unit_cell()) {
        trajan::log::warn("wrap: frame has no unit cell, skipping");
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
  } else {
    trajan::log::debug("wrap: registering molecular PBC wrap");
    pipeline.add_transform("wrap", [&traj](trajan::core::Frame &frame) mutable {
      if (!frame.has_unit_cell()) {
        trajan::log::warn("wrap: frame has no unit cell, skipping");
        return;
      }
      const auto &uc = frame.unit_cell().value();
      const auto molecules = traj.get_molecules();
      for (const auto &mol : molecules) {
        occ::Vec3 com = mol.center_of_mass();
        occ::Vec3 frac_com = uc.to_fractional(com);
        occ::Vec3 image = frac_com.array().floor();
        occ::Vec3 cart_shift = uc.to_cartesian(image);
        for (const auto &atom : mol.atoms()) {
          frame.atoms()[atom.index].x -= cart_shift[0];
          frame.atoms()[atom.index].y -= cart_shift[1];
          frame.atoms()[atom.index].z -= cart_shift[2];
        }
      }
    });
  }
}

CLI::App *add_wrap_subcommand(CLI::App &app, Trajectory &traj,
                              Pipeline &pipeline) {
  CLI::App *wrap = app.add_subcommand(
      "wrap", "Wrap atoms/molecules into the primary (0,0,0) unit cell image. "
              "Requires load subcommand.");

  auto opts = std::make_shared<WrapOpts>();
  wrap->add_flag("--atomic", opts->atomic,
                 "Wrap individual atoms rather than whole molecules. By "
                 "default molecules are wrapped using their center of mass.");

  wrap->callback([opts, &traj, &pipeline]() {
    trajan::log::set_subcommand_log_pattern("wrap");
    trajan::log::debug("Beginning wrap subcommand");
    register_wrap_transform(*opts, traj, pipeline);
  });

  return wrap;
}

} // namespace trajan::main
