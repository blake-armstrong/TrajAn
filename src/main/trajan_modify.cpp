#include <ankerl/unordered_dense.h>
#include <cmath>
#include <stdexcept>
#include <trajan/core/log.h>
#include <trajan/core/neigh.h>
#include <trajan/core/units.h>
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

  if (!opts.raw_rotate.empty()) {
    double angle_deg = std::stod(opts.raw_rotate[0]);
    char axis_char =
        static_cast<char>(std::tolower(static_cast<unsigned char>(
            opts.raw_rotate[1][0])));
    if (axis_char != 'x' && axis_char != 'y' && axis_char != 'z') {
      throw std::invalid_argument(
          fmt::format("--rotate axis must be x, y, or z, got '{}'",
                      opts.raw_rotate[1]));
    }
    double angle_rad = trajan::units::radians(angle_deg);
    double c = std::cos(angle_rad), s = std::sin(angle_rad);
    double ox = opts.rotate_origin[0];
    double oy = opts.rotate_origin[1];
    double oz = opts.rotate_origin[2];

    // Rotate about an arbitrary point: translate to origin, rotate, translate back.
    auto rotate_pos = [axis_char, c, s, ox, oy, oz](double &x, double &y,
                                                      double &z) {
      x -= ox; y -= oy; z -= oz;
      double nx, ny, nz;
      switch (axis_char) {
      case 'x':
        nx = x; ny = c * y - s * z; nz = s * y + c * z;
        break;
      case 'y':
        nx = c * x + s * z; ny = y; nz = -s * x + c * z;
        break;
      default: // z
        nx = c * x - s * y; ny = s * x + c * y; nz = z;
        break;
      }
      x = nx + ox; y = ny + oy; z = nz + oz;
    };

    if (opts.raw_sel.empty()) {
      trajan::log::debug(
          "modify: registering rotate {:.4f} deg about {} for all atoms",
          angle_deg, axis_char);
      pipeline.add_transform(
          [rotate_pos](trajan::core::Frame &frame) mutable {
            for (auto &atom : frame.atoms())
              rotate_pos(atom.x, atom.y, atom.z);
          });
    } else {
      auto parsed_sel = io::selection_expr_validator(opts.raw_sel);
      io::MolOrigin mol_origin = opts.mol_origin;
      trajan::log::debug(
          "modify: registering rotate {:.4f} deg about {} for selection '{}'",
          angle_deg, axis_char, opts.raw_sel);

      bool first_call = true;
      std::vector<core::EntityVariant> entities;
      pipeline.add_transform([rotate_pos, parsed_sel, mol_origin, &traj,
                               first_call,
                               entities](trajan::core::Frame &frame) mutable {
        if (first_call || traj.topology_has_changed()) {
          entities = traj.get_entities(parsed_sel, mol_origin);
          first_call = false;
        } else {
          traj.update_entities(entities);
        }
        ankerl::unordered_dense::set<size_t> indices;
        indices.reserve(entities.size());
        for (const auto &ev : entities) {
          std::visit(
              [&](const auto &e) {
                using T = std::decay_t<decltype(e)>;
                if constexpr (std::is_same_v<T, core::Atom>) {
                  indices.insert(static_cast<size_t>(e.index));
                } else if constexpr (std::is_same_v<T, core::Molecule>) {
                  for (const auto &a : e.atoms())
                    indices.insert(static_cast<size_t>(a.index));
                }
              },
              ev);
        }
        for (auto &atom : frame.atoms()) {
          if (indices.count(static_cast<size_t>(atom.index)))
            rotate_pos(atom.x, atom.y, atom.z);
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

  // Helper: build an index set from entities (reused by all three setters).
  auto make_index_set =
      [](const std::vector<core::EntityVariant> &entities) {
        ankerl::unordered_dense::set<size_t> indices;
        indices.reserve(entities.size());
        for (const auto &ev : entities) {
          std::visit(
              [&](const auto &e) {
                using T = std::decay_t<decltype(e)>;
                if constexpr (std::is_same_v<T, core::Atom>) {
                  indices.insert(static_cast<size_t>(e.index));
                } else if constexpr (std::is_same_v<T, core::Molecule>) {
                  for (const auto &a : e.atoms())
                    indices.insert(static_cast<size_t>(a.index));
                }
              },
              ev);
        }
        return indices;
      };

  if (opts.set_charge.has_value()) {
    const double val = opts.set_charge.value();
    if (opts.raw_sel.empty()) {
      trajan::log::debug("modify: set charge = {} for all atoms", val);
      pipeline.add_transform([val](trajan::core::Frame &frame) {
        for (auto &atom : frame.atoms())
          atom.charge = val;
      });
    } else {
      auto parsed_sel = io::selection_expr_validator(opts.raw_sel);
      io::MolOrigin mol_origin = opts.mol_origin;
      trajan::log::debug("modify: set charge = {} for selection '{}'", val,
                         opts.raw_sel);
      bool first_call = true;
      std::vector<core::EntityVariant> entities;
      pipeline.add_transform(
          [val, parsed_sel, mol_origin, make_index_set, &traj, first_call,
           entities](trajan::core::Frame &frame) mutable {
            if (first_call || traj.topology_has_changed()) {
              entities = traj.get_entities(parsed_sel, mol_origin);
              first_call = false;
            } else {
              traj.update_entities(entities);
            }
            auto indices = make_index_set(entities);
            for (auto &atom : frame.atoms())
              if (indices.count(static_cast<size_t>(atom.index)))
                atom.charge = val;
          });
    }
  }

  if (opts.set_type.has_value()) {
    const std::string val = opts.set_type.value();
    if (opts.raw_sel.empty()) {
      trajan::log::debug("modify: set type = '{}' for all atoms", val);
      pipeline.add_transform([val](trajan::core::Frame &frame) {
        for (auto &atom : frame.atoms())
          atom.type = val;
      });
    } else {
      auto parsed_sel = io::selection_expr_validator(opts.raw_sel);
      io::MolOrigin mol_origin = opts.mol_origin;
      trajan::log::debug("modify: set type = '{}' for selection '{}'", val,
                         opts.raw_sel);
      bool first_call = true;
      std::vector<core::EntityVariant> entities;
      pipeline.add_transform(
          [val, parsed_sel, mol_origin, make_index_set, &traj, first_call,
           entities](trajan::core::Frame &frame) mutable {
            if (first_call || traj.topology_has_changed()) {
              entities = traj.get_entities(parsed_sel, mol_origin);
              first_call = false;
            } else {
              traj.update_entities(entities);
            }
            auto indices = make_index_set(entities);
            for (auto &atom : frame.atoms())
              if (indices.count(static_cast<size_t>(atom.index)))
                atom.type = val;
          });
    }
  }

  if (opts.set_mol_type.has_value()) {
    const std::string val = opts.set_mol_type.value();
    if (opts.raw_sel.empty()) {
      trajan::log::debug("modify: set mol-type = '{}' for all atoms", val);
      pipeline.add_transform([val](trajan::core::Frame &frame) {
        for (auto &atom : frame.atoms())
          atom.molecule_type = val;
      });
    } else {
      auto parsed_sel = io::selection_expr_validator(opts.raw_sel);
      io::MolOrigin mol_origin = opts.mol_origin;
      trajan::log::debug("modify: set mol-type = '{}' for selection '{}'", val,
                         opts.raw_sel);
      bool first_call = true;
      std::vector<core::EntityVariant> entities;
      pipeline.add_transform(
          [val, parsed_sel, mol_origin, make_index_set, &traj, first_call,
           entities](trajan::core::Frame &frame) mutable {
            if (first_call || traj.topology_has_changed()) {
              entities = traj.get_entities(parsed_sel, mol_origin);
              first_call = false;
            } else {
              traj.update_entities(entities);
            }
            auto indices = make_index_set(entities);
            for (auto &atom : frame.atoms())
              if (indices.count(static_cast<size_t>(atom.index)))
                atom.molecule_type = val;
          });
    }
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

  auto *rotate_opt =
      modify->add_option("--rotate", opts->raw_rotate,
                         "Rotate atom positions by <degrees> about <axis> "
                         "(axis: x, y, or z). Affects all atoms unless --sel "
                         "is given.\nExample: --rotate 45.0 z")
          ->expected(2);

  modify->add_option("--rotate-origin", opts->rotate_origin,
                     "Point to rotate about (default: 0 0 0).\n"
                     "Example: --rotate-origin 5.0 5.0 0.0")
      ->expected(3)
      ->capture_default_str();

  modify->add_flag("--wrap", opts->wrap,
                   "Wrap atom positions into the primary unit cell (requires "
                   "unit cell in trajectory)");

  modify->add_option("--set-charge", opts->set_charge,
                     "Set the partial charge of selected atoms to this value. "
                     "Affects all atoms unless --sel is given.");

  modify->add_option("--set-type", opts->set_type,
                     "Set the atom type label of selected atoms (e.g. CA, OW). "
                     "Affects all atoms unless --sel is given.");

  modify->add_option("--set-mol-type", opts->set_mol_type,
                     "Set the molecule/residue type of selected atoms (e.g. WAT, ASP). "
                     "Affects all atoms unless --sel is given.");

  modify->add_option(
      "--sel,-s", opts->raw_sel,
      "Selection to restrict --translate and --rotate (prefix: i=atom indices,\n"
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

  modify->callback([opts, trans_opt, rotate_opt, &traj, &pipeline]() {
    trajan::log::set_subcommand_log_pattern("modify");
    trajan::log::debug("Beginning modify subcommand");

    if (trans_opt->count() == 0) {
      opts->translate = {0.0, 0.0, 0.0};
    }
    if (rotate_opt->count() == 0) {
      opts->raw_rotate.clear();
    }

    register_modify_transforms(*opts, traj, pipeline);
    trajan::log::debug("modify subcommand registered {} transform(s)",
                       pipeline.size());
  });

  return modify;
}

} // namespace trajan::main
