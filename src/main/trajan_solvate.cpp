#include <trajan/core/frame.h>
#include <trajan/core/log.h>
#include <trajan/core/neigh.h>
#include <trajan/core/trajectory.h>
#include <trajan/io/file_handler.h>
#include <trajan/main/trajan_solvate.h>

#include <occ/core/linear_algebra.h>
#include <occ/crystal/unitcell.h>

#include <cmath>
#include <limits>
#include <numbers>
#include <random>
#include <stdexcept>

namespace trajan::main {

static occ::crystal::UnitCell
parse_box(const std::vector<double> &box) {
  switch (box.size()) {
  case 1:
    return occ::crystal::orthorhombic_cell(box[0], box[0], box[0]);
  case 3:
    return occ::crystal::orthorhombic_cell(box[0], box[1], box[2]);
  case 6: {
    const double deg2rad = std::numbers::pi / 180.0;
    return occ::crystal::triclinic_cell(box[0], box[1], box[2],
                                        box[3] * deg2rad, box[4] * deg2rad,
                                        box[5] * deg2rad);
  }
  default:
    throw std::runtime_error(
        "solvate: --box requires 1 (cubic), 3 (orthorhombic), or 6 "
        "(triclinic a b c α β γ) values");
  }
}

static occ::Vec3 cartesian_bounding_box(const occ::crystal::UnitCell &uc) {
  const occ::Mat3 &D = uc.direct();
  occ::Mat3N corners(3, 8);
  corners.col(0) = occ::Vec3::Zero();
  corners.col(1) = D.col(0);
  corners.col(2) = D.col(1);
  corners.col(3) = D.col(2);
  corners.col(4) = D.col(0) + D.col(1);
  corners.col(5) = D.col(0) + D.col(2);
  corners.col(6) = D.col(1) + D.col(2);
  corners.col(7) = D.col(0) + D.col(1) + D.col(2);
  return corners.rowwise().maxCoeff() - corners.rowwise().minCoeff();
}

static occ::Mat3 random_rotation(std::mt19937 &rng) {
  std::normal_distribution<double> normal(0.0, 1.0);
  double w = normal(rng), x = normal(rng), y = normal(rng), z = normal(rng);
  const double n = std::sqrt(w*w + x*x + y*y + z*z);
  w /= n; x /= n; y /= n; z /= n;
  occ::Mat3 R;
  R << 1 - 2*(y*y + z*z),     2*(x*y - w*z),     2*(x*z + w*y),
           2*(x*y + w*z), 1 - 2*(x*x + z*z),     2*(y*z - w*x),
           2*(x*z - w*y),     2*(y*z + w*x), 1 - 2*(x*x + y*y);
  return R;
}

void run_solvate_subcommand(const SolvateOpts &opts, Trajectory &traj,
                            Pipeline &pipeline) {
  auto tmpl_handler = trajan::io::read_input_file(opts.molecule_file);
  tmpl_handler->initialise(trajan::io::FileHandler::Mode::Read);
  trajan::core::Frame tmpl_frame;
  if (!tmpl_handler->read_frame(tmpl_frame))
    throw std::runtime_error("solvate: could not read template molecule file");
  tmpl_handler->finalise();

  const auto tmpl_atoms = tmpl_frame.atoms();
  if (tmpl_atoms.empty())
    throw std::runtime_error("solvate: template molecule has no atoms");

  const int n_tmpl = static_cast<int>(tmpl_atoms.size());
  occ::Mat3N tmpl_pos(3, n_tmpl);
  for (int i = 0; i < n_tmpl; ++i) {
    tmpl_pos(0, i) = tmpl_atoms[i].x;
    tmpl_pos(1, i) = tmpl_atoms[i].y;
    tmpl_pos(2, i) = tmpl_atoms[i].z;
  }
  tmpl_pos.colwise() -= tmpl_pos.rowwise().mean().eval();

  double molar_mass = 0.0;
  for (const auto &a : tmpl_atoms)
    molar_mass += a.element.mass();

  // Extract unique bonds from the template graph (vertex descriptor == local atom index).
  struct TmplBond { int i, j; double length; };
  std::vector<TmplBond> tmpl_bonds;
  {
    const auto &g = tmpl_frame.get_atom_graph();
    for (const auto &[v, neighbors] : g.adjacency_list()) {
      for (const auto &[nb, ed] : neighbors) {
        if (v < nb)
          tmpl_bonds.push_back({static_cast<int>(v), static_cast<int>(nb),
                                 g.edge(ed).bond_length});
      }
    }
  }
  if (!tmpl_bonds.empty())
    trajan::log::debug("solvate: template has {} bonds", tmpl_bonds.size());

  pipeline.add_transform("solvate", [opts, tmpl_atoms, tmpl_pos, molar_mass,
                          n_tmpl, tmpl_bonds](trajan::core::Frame &frame) mutable {
    const occ::Vec3 origin(opts.origin[0], opts.origin[1], opts.origin[2]);
    const occ::crystal::UnitCell box_uc = parse_box(opts.box);
    const occ::Vec3 bbox = cartesian_bounding_box(box_uc);
    const double volume = box_uc.volume();

    trajan::core::CellList cl(opts.min_dist);
    cl.initialize_region(origin, bbox, box_uc);

    std::mt19937 rng(opts.seed);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    int n_target = std::numeric_limits<int>::max();
    if (opts.density.has_value()) {
      n_target = static_cast<int>(
          std::round(opts.density.value() * volume / (molar_mass * 1.6605388)));
      trajan::log::info(
          "solvate: molar mass = {:.4f} g/mol, fill volume = {:.2f} Å³",
          molar_mass, volume);
      trajan::log::info("solvate: targeting {} molecules for {:.4f} g/cm³",
                        n_target, opts.density.value());
    } else if (opts.n_molecules.has_value()) {
      n_target = opts.n_molecules.value();
    }

    const int max_total =
        (opts.n_molecules.has_value() || opts.density.has_value())
            ? n_target * opts.max_attempts
            : opts.max_attempts;

    int global_atom_idx = static_cast<int>(frame.atoms().size());
    int mol_idx = 0;
    for (const auto &a : frame.atoms())
      mol_idx = std::max(mol_idx, a.molecule_index + 1);

    std::vector<trajan::core::EnhancedAtom> placed_atoms;
    std::vector<int> copy_starts;
    int placed = 0;
    int attempts = 0;

    for (; attempts < max_total && placed < n_target; ++attempts) {
      const occ::Vec3 frac(uniform(rng), uniform(rng), uniform(rng));
      const occ::Vec3 pos = origin + box_uc.direct() * frac;

      const occ::Mat3 R = random_rotation(rng);
      const occ::Mat3N candidate = (R * tmpl_pos).colwise() + pos;

      bool clash = false;
      for (int i = 0; i < n_tmpl && !clash; ++i)
        clash = cl.has_clash(candidate.col(i));
      if (clash)
        continue;

      const int copy_start = global_atom_idx;
      for (int i = 0; i < n_tmpl; ++i) {
        cl.insert_point(candidate.col(i), global_atom_idx);
        trajan::core::EnhancedAtom atom = tmpl_atoms[i];
        atom.x = candidate(0, i);
        atom.y = candidate(1, i);
        atom.z = candidate(2, i);
        atom.index = global_atom_idx++;
        atom.molecule_index = mol_idx;
        placed_atoms.push_back(atom);
      }
      copy_starts.push_back(copy_start);
      ++mol_idx;
      ++placed;
    }

    if (placed < n_target)
      trajan::log::warn("solvate: placed {}/{} molecules after {} attempts",
                        placed, n_target, attempts);
    else
      trajan::log::info("solvate: placed {} molecules in {} attempts", placed,
                        attempts);

    auto all_atoms = frame.atoms();
    const size_t slab_count = all_atoms.size();
    all_atoms.insert(all_atoms.end(), placed_atoms.begin(), placed_atoms.end());
    frame.set_atoms(all_atoms);
    trajan::log::info("solvate: {} total atoms ({} original + {} placed)",
                      all_atoms.size(), slab_count, placed_atoms.size());

    if (!tmpl_bonds.empty()) {
      auto graph = frame.get_atom_graph();
      for (int copy_start : copy_starts) {
        for (int i = 0; i < n_tmpl; ++i)
          graph.add_vertex(trajan::core::AtomVertex{copy_start + i});
        for (const auto &b : tmpl_bonds) {
          trajan::core::BondEdge edge;
          edge.indices = {copy_start + b.i, copy_start + b.j};
          edge.bond_length = b.length;
          graph.add_edge(copy_start + b.i, copy_start + b.j, edge, true);
        }
      }
      frame.set_atom_graph(graph);
    }
  });
}

CLI::App *add_solvate_subcommand(CLI::App &app, Trajectory &traj,
                                 Pipeline &pipeline) {
  CLI::App *solv = app.add_subcommand(
      "solvate",
      "Fill a region of space with copies of a template molecule using "
      "random sequential adsorption. Requires load subcommand.");

  auto opts = std::make_shared<SolvateOpts>();

  solv->add_option("-m,--molecule", opts->molecule_file,
                   "Template molecule file (PDB, XYZ, etc.)")
      ->required();
  solv->add_option("--origin", opts->origin,
                   "Cartesian origin of the fill region in Å (x y z)")
      ->expected(3)
      ->default_str("0 0 0");
  solv->add_option("--box", opts->box,
                   "Fill region dimensions: 1 value (cubic a), "
                   "3 values (orthorhombic a b c), or "
                   "6 values (triclinic a b c α β γ with angles in degrees)")
      ->expected(1, 6)
      ->required();
  auto *n_opt = solv->add_option("-n,--n-molecules", opts->n_molecules,
                                 "Number of molecules to place.");
  auto *density_opt = solv->add_option(
      "--density", opts->density,
      "Target mass density (g/cm³). Computes the required number of molecules "
      "from the fill volume and template molar mass.");
  n_opt->excludes(density_opt);
  density_opt->excludes(n_opt);
  solv->add_option("--min-dist", opts->min_dist,
                   "Minimum allowed distance between any two atoms (Å)")
      ->default_val(2.5);
  solv->add_option("--max-attempts", opts->max_attempts,
                   "Max placement attempts per molecule (with -n or --density),"
                   " or total attempts otherwise.")
      ->default_val(200);
  solv->add_option("--seed", opts->seed, "Random seed")->default_val(42);

  solv->callback([opts, &traj, &pipeline]() {
    trajan::log::set_subcommand_log_pattern("solvate");
    trajan::log::debug("Beginning solvate subcommand");
    run_solvate_subcommand(*opts, traj, pipeline);
  });

  return solv;
}

} // namespace trajan::main
