#include <trajan/core/frame.h>
#include <trajan/core/log.h>
#include <trajan/core/trajectory.h>
#include <trajan/io/file_handler.h>
#include <trajan/main/trajan_surface.h>

#include <occ/crystal/asymmetric_unit.h>
#include <occ/crystal/crystal.h>
#include <occ/crystal/hkl.h>
#include <occ/crystal/spacegroup.h>
#include <occ/crystal/surface.h>
#include <occ/crystal/unitcell.h>

#include <algorithm>
#include <stdexcept>

namespace trajan::main {

using occ::crystal::AsymmetricUnit;
using occ::crystal::Crystal;
using occ::crystal::HKL;
using occ::crystal::SpaceGroup;
using occ::crystal::Surface;
using occ::crystal::UnitCell;
using occ::crystal::triclinic_cell;

// Build an occ Crystal in P1 symmetry from a trajan Frame.
// The Crystal only needs its unit cell for Surface construction, so
// set_connectivity_criteria(false) avoids unnecessary bond detection.
static Crystal build_p1_crystal(const trajan::core::Frame &frame) {
  const auto &atoms = frame.atoms();
  const UnitCell &uc = frame.unit_cell().value();
  const size_t n = atoms.size();

  occ::Mat3N frac_pos(3, n);
  occ::IVec atomic_numbers(n);
  std::vector<std::string> labels(n);

  for (size_t i = 0; i < n; ++i) {
    occ::Vec3 cart(atoms[i].x, atoms[i].y, atoms[i].z);
    frac_pos.col(i) = uc.to_fractional(cart);
    atomic_numbers(i) = atoms[i].element.atomic_number();
    labels[i] = atoms[i].type;
  }

  AsymmetricUnit asym(frac_pos, atomic_numbers, labels);
  SpaceGroup sg; // P1 — all atoms in the asymmetric unit
  Crystal crystal(asym, sg, uc);
  crystal.set_connectivity_criteria(false);
  return crystal;
}

// Derive a UnitCell from the surface basis matrix.
static UnitCell surface_unit_cell(const Surface &surf, double depth_scale) {
  const occ::Mat3 basis = surf.basis_matrix(depth_scale);
  const occ::Vec3 a = basis.col(0);
  const occ::Vec3 b = basis.col(1);
  const occ::Vec3 c = basis.col(2);
  const double la = a.norm(), lb = b.norm(), lc = c.norm();
  const double cos_alpha = b.normalized().dot(c.normalized());
  const double cos_beta  = a.normalized().dot(c.normalized());
  const double cos_gamma = a.normalized().dot(b.normalized());
  const double alpha = std::acos(std::clamp(cos_alpha, -1.0, 1.0));
  const double beta  = std::acos(std::clamp(cos_beta,  -1.0, 1.0));
  const double gamma = std::acos(std::clamp(cos_gamma, -1.0, 1.0));
  return triclinic_cell(la, lb, lc, alpha, beta, gamma);
}

// Build a trajan Frame from molecules returned by find_molecule_cell_translations
// (molecular mode — molecules kept intact).
static trajan::core::Frame
build_mol_frame(const std::vector<occ::core::Molecule> &slab_mols,
                const std::vector<trajan::core::EnhancedMolecule> &orig_mols,
                const UnitCell &slab_uc) {
  trajan::core::Frame frame;
  frame.set_unit_cell(slab_uc);

  std::vector<trajan::core::EnhancedAtom> out_atoms;
  int mol_idx = 0;
  for (const auto &mol : slab_mols) {
    const int uc_mol_idx = mol.unit_cell_molecule_idx();
    const auto &orig_mol = orig_mols[uc_mol_idx];
    const auto &pos = mol.positions(); // 3×N Angstroms

    for (size_t ai = 0; ai < orig_mol.enhanced_atoms.size(); ++ai) {
      trajan::core::EnhancedAtom atom = orig_mol.enhanced_atoms[ai];
      atom.x = pos(0, ai);
      atom.y = pos(1, ai);
      atom.z = pos(2, ai);
      atom.index = static_cast<int>(out_atoms.size());
      atom.molecule_index = mol_idx;
      out_atoms.push_back(atom);
    }
    ++mol_idx;
  }

  frame.set_atoms(out_atoms);
  return frame;
}

// Build a trajan Frame from single-atom molecules (atomic mode).
static trajan::core::Frame
build_atomic_frame(const std::vector<occ::core::Molecule> &slab_mols,
                   const std::vector<trajan::core::EnhancedAtom> &orig_atoms,
                   const UnitCell &slab_uc) {
  trajan::core::Frame frame;
  frame.set_unit_cell(slab_uc);

  std::vector<trajan::core::EnhancedAtom> out_atoms;
  for (const auto &mol : slab_mols) {
    const int orig_idx = mol.unit_cell_molecule_idx();
    trajan::core::EnhancedAtom atom = orig_atoms[orig_idx];
    const auto &pos = mol.positions();
    atom.x = pos(0, 0);
    atom.y = pos(1, 0);
    atom.z = pos(2, 0);
    atom.index = static_cast<int>(out_atoms.size());
    out_atoms.push_back(atom);
  }

  frame.set_atoms(out_atoms);
  return frame;
}

void run_surface_subcommand(const SurfaceOpts &opts, Trajectory &traj) {
  if (!traj.next_frame())
    throw std::runtime_error("surface: no frames available from input");

  trajan::core::Frame &frame = traj.frame();

  if (!frame.has_unit_cell())
    throw std::runtime_error(
        "surface: input has no unit cell — a CRYST1 record is required");

  const UnitCell &uc = frame.unit_cell().value();
  if (opts.hkl.size() != 3)
    throw std::runtime_error("surface: hkl must have exactly 3 indices");

  const HKL hkl{opts.hkl[0], opts.hkl[1], opts.hkl[2]};

  const Crystal crystal = build_p1_crystal(frame);
  Surface surface(hkl, crystal);

  trajan::log::info("Surface ({} {} {}): d-spacing = {:.3f} Å, area = {:.3f} Å²",
                    opts.hkl[0], opts.hkl[1], opts.hkl[2],
                    surface.d(), surface.area());

  const UnitCell slab_uc = surface_unit_cell(surface, opts.depth);

  // Build the list of unit-cell entities (molecules or atoms) to cut.
  std::vector<occ::core::Molecule> uc_mols;
  std::vector<trajan::core::EnhancedMolecule> enhanced_mols;
  std::vector<trajan::core::EnhancedAtom> frame_atoms;

  if (opts.atomic) {
    frame_atoms = frame.atoms();
    for (size_t i = 0; i < frame_atoms.size(); ++i) {
      occ::IVec nums(1);
      nums(0) = frame_atoms[i].element.atomic_number();
      occ::Mat3N p(3, 1);
      p(0, 0) = frame_atoms[i].x;
      p(1, 0) = frame_atoms[i].y;
      p(2, 0) = frame_atoms[i].z;
      occ::core::Molecule mol(nums, p);
      mol.set_unit_cell_molecule_idx(i);
      uc_mols.push_back(std::move(mol));
    }
    trajan::log::info("Atomic mode: {} unit-cell atoms", uc_mols.size());
  } else {
    enhanced_mols = traj.get_molecules();
    uc_mols.reserve(enhanced_mols.size());
    for (size_t i = 0; i < enhanced_mols.size(); ++i) {
      enhanced_mols[i].sync_occ_positions();
      // Slice to base occ::core::Molecule (positions already synced).
      occ::core::Molecule base_mol =
          static_cast<const occ::core::Molecule &>(enhanced_mols[i]);
      base_mol.set_unit_cell_molecule_idx(i);
      uc_mols.push_back(std::move(base_mol));
    }
    trajan::log::info("Molecular mode: {} unit-cell molecules", uc_mols.size());
  }

  // Centroids used to compute symmetry-unique cuts.
  occ::Mat3N centroids(3, uc_mols.size());
  for (size_t i = 0; i < uc_mols.size(); ++i)
    centroids.col(i) = uc_mols[i].centroid();

  // Determine which cut offsets to use.
  std::vector<double> cuts;
  if (opts.shift.has_value()) {
    cuts.push_back(opts.shift.value());
  } else {
    cuts = surface.possible_cuts(centroids);
    trajan::log::info("Found {} possible cut position(s)", cuts.size());
  }

  auto out_handler = trajan::io::write_output_file(opts.outfile);

  int frame_count = 0;
  for (double cut : cuts) {
    trajan::log::info("Cut shift = {:.6f}", cut);
    auto slab_mols =
        surface.find_molecule_cell_translations(uc_mols, opts.depth, cut);

    trajan::core::Frame slab_frame =
        opts.atomic ? build_atomic_frame(slab_mols, frame_atoms, slab_uc)
                    : build_mol_frame(slab_mols, enhanced_mols, slab_uc);

    trajan::log::info("  frame {}: {} atoms", frame_count + 1,
                      slab_frame.num_atoms());
    out_handler->write_frame(slab_frame);
    ++frame_count;
  }

  out_handler->finalise();
  trajan::log::info("Wrote {} surface frame(s) to {}", frame_count,
                    opts.outfile.string());
}

CLI::App *add_surface_subcommand(CLI::App &app, Trajectory &traj,
                                 Pipeline &pipeline) {
  CLI::App *surf = app.add_subcommand(
      "surface",
      "Cut a surface slab from a bulk crystal using Miller indices (h k l). "
      "Without --shift, all symmetry-unique cuts are generated as separate "
      "frames in a PDB trajectory. Requires load subcommand.");

  auto opts = std::make_shared<SurfaceOpts>();
  surf->add_option("hkl", opts->hkl,
                   "Miller indices h k l (three integers, e.g. 1 0 0)")
      ->required()
      ->expected(3);
  surf->add_option("outfile", opts->outfile, "Output PDB trajectory file")
      ->required();
  surf->add_option("--shift", opts->shift,
                   "Fractional shift along the surface normal (0–1). "
                   "If omitted, all possible cuts are generated.");
  surf->add_option("--depth", opts->depth,
                   "Slab depth as a multiple of the interplanar d-spacing")
      ->default_val(1.0);
  surf->add_flag("--atomic", opts->atomic,
                 "Cut at the atom level (molecules may be broken). "
                 "By default whole molecules are kept intact.");

  surf->callback([opts, &traj, &pipeline]() {
    trajan::log::set_subcommand_log_pattern("surface");
    trajan::log::debug("Beginning surface subcommand");
    run_surface_subcommand(*opts, traj);
  });

  return surf;
}

} // namespace trajan::main
