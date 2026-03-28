#include <ankerl/unordered_dense.h>
#include <trajan/core/log.h>
#include <trajan/core/neigh.h>
#include <trajan/io/selection.h>
#include <trajan/main/trajan_modify.h>

namespace trajan::main {

void register_modify_transforms(const ModifyOpts &opts, Trajectory &traj,
                                Pipeline &pipeline) {
  const double dx = opts.translate[0];
  const double dy = opts.translate[1];
  const double dz = opts.translate[2];

  if (dx != 0.0 || dy != 0.0 || dz != 0.0) {
    if (opts.raw_sel.empty()) {
      trajan::log::debug("modify: registering translate ({}, {}, {}) for all "
                         "atoms",
                         dx, dy, dz);
      pipeline.add_transform([dx, dy, dz](trajan::core::Frame &frame) {
        for (auto &atom : frame.atoms()) {
          atom.x += dx;
          atom.y += dy;
          atom.z += dz;
        }
      });
    } else {
      auto parsed_sel = io::selection_expr_validator(opts.raw_sel);
      io::MolOrigin mol_origin = opts.mol_origin;
      trajan::log::debug(
          "modify: registering translate ({}, {}, {}) for selection '{}'", dx,
          dy, dz, opts.raw_sel);

      bool first_call = true;
      std::vector<core::EntityVariant> entities;
      pipeline.add_transform([dx, dy, dz, parsed_sel, mol_origin, &traj,
                               first_call,
                               entities](trajan::core::Frame &frame) mutable {
        if (first_call || traj.topology_has_changed()) {
          entities = traj.get_entities(parsed_sel, mol_origin);
          first_call = false;
        } else {
          traj.update_entities(entities);
        }

        // Collect atom indices to translate. For molecule entities, include
        // all atoms belonging to the molecule.
        ankerl::unordered_dense::set<size_t> indices;
        indices.reserve(entities.size());
        for (const auto &ev : entities) {
          std::visit(
              [&](const auto &e) {
                using T = std::decay_t<decltype(e)>;
                if constexpr (std::is_same_v<T, core::Atom>) {
                  indices.insert(static_cast<size_t>(e.index));
                } else if constexpr (std::is_same_v<T, core::Molecule>) {
                  for (const auto &a : e.atoms()) {
                    indices.insert(static_cast<size_t>(a.index));
                  }
                }
              },
              ev);
        }

        for (auto &atom : frame.atoms()) {
          if (indices.count(static_cast<size_t>(atom.index))) {
            atom.x += dx;
            atom.y += dy;
            atom.z += dz;
          }
        }
      });
    }
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
                         "Translate atom positions by dx dy dz (Angstroms). "
                         "Affects all atoms unless --sel is given.")
          ->expected(3);

  modify->add_flag("--wrap", opts->wrap,
                   "Wrap atom positions into the primary unit cell (requires "
                   "unit cell in trajectory)");

  modify->add_option(
      "--sel,-s", opts->raw_sel,
      "Selection to restrict --translate (prefix: i=atom indices, "
      "a=atom types, j=molecule indices, m=molecule types).\n"
      "Supports boolean expressions and position predicates, e.g.:\n"
      "  (aO) and (z > 40 and z < 80)\n"
      "  not (aH)\n"
      "  mM1 and z > 10");

  modify->add_option(
      "--mol-origin", opts->mol_origin,
      "Reference position used for position predicates on molecule selections.\n"
      "  com      - center of mass (default)\n"
      "  centroid - unweighted centroid")
      ->transform(CLI::CheckedTransformer(
          std::map<std::string, io::MolOrigin>{
              {"com", io::MolOrigin::CenterOfMass},
              {"centroid", io::MolOrigin::Centroid}}))
      ->capture_default_str();

  modify->callback([opts, trans_opt, &traj, &pipeline]() {
    trajan::log::set_subcommand_log_pattern("modify");
    trajan::log::debug("Beginning modify subcommand");

    if (trans_opt->count() == 0) {
      opts->translate = {0.0, 0.0, 0.0};
    }

    register_modify_transforms(*opts, traj, pipeline);
    trajan::log::debug("modify subcommand registered {} transform(s)",
                       pipeline.size());
  });

  return modify;
}

} // namespace trajan::main
