#include <fstream>
#include <cmath>
#include <stdexcept>
#include <trajan/core/log.h>
#include <trajan/core/topology.h>
#include <trajan/main/trajan_dipole.h>
#include <occ/core/constants.h>
#include <occ/crystal/unitcell.h>

namespace trajan::main {

using occ::crystal::UnitCell;

// ── B-spline helpers ───────────────────────────────────────────────────────

static double bspline(int p, double x) {
  if (x < 0.0 || x >= static_cast<double>(p))
    return 0.0;
  if (p == 1)
    return 1.0;
  return (x / (p - 1)) * bspline(p - 1, x) +
         ((p - x) / (p - 1)) * bspline(p - 1, x - 1.0);
}

static std::vector<double> bspline_weights(double w, int p) {
  const int i0 = static_cast<int>(std::floor(w));
  const double frac = w - i0;
  std::vector<double> wt(p);
  for (int k = 0; k < p; ++k)
    wt[k] = bspline(p, frac + k);
  return wt;
}

static size_t fft_friendly(size_t n) {
  if (n < 1) n = 1;
  while (true) {
    size_t t = n;
    while (t % 2 == 0) t /= 2;
    while (t % 3 == 0) t /= 3;
    while (t % 5 == 0) t /= 5;
    if (t == 1) return n;
    ++n;
  }
}

// ── DX output ─────────────────────────────────────────────────────────────

static void write_dipole_dx(const std::string &path,
                            const std::vector<double> &data,
                            size_t Na, size_t Nb, size_t Nc,
                            const UnitCell &uc,
                            const std::string &comment) {
  std::ofstream f(path);
  if (!f)
    throw std::runtime_error("dipole: cannot open output file: " + path);

  const occ::Vec3 da = uc.a_vector() / static_cast<double>(Na);
  const occ::Vec3 db = uc.b_vector() / static_cast<double>(Nb);
  const occ::Vec3 dc = uc.c_vector() / static_cast<double>(Nc);

  f << "# " << comment << "\n";
  f << "object 1 class gridpositions counts " << Na << " " << Nb << " " << Nc << "\n";
  f << "origin 0.0 0.0 0.0\n";
  f << "delta " << da[0] << " " << da[1] << " " << da[2] << "\n";
  f << "delta " << db[0] << " " << db[1] << " " << db[2] << "\n";
  f << "delta " << dc[0] << " " << dc[1] << " " << dc[2] << "\n";
  f << "object 2 class gridconnections counts " << Na << " " << Nb << " " << Nc << "\n";
  f << "object 3 class array type double rank 0 items " << (Na * Nb * Nc)
    << " data follows\n";

  size_t col = 0;
  for (const double v : data) {
    f << v;
    if (++col % 3 == 0) f << "\n";
    else f << "\t";
  }
  if (col % 3 != 0) f << "\n";

  f << "object \"" << comment << "\" class field\n"
    << "component \"positions\" value 1\n"
    << "component \"connections\" value 2\n"
    << "component \"data\" value 3\n"
    << "attribute \"dep\" string \"positions\"\n";
}

// ── Polarization grid ──────────────────────────────────────────────────────

struct DipoleGrid {
  size_t Na, Nb, Nc;
  int order;

  std::vector<double> Px, Py, Pz;  // accumulated dipole components (e·Å)
  size_t frame_count{0};

  DipoleGrid(size_t Na, size_t Nb, size_t Nc, int order)
      : Na(Na), Nb(Nb), Nc(Nc), order(order),
        Px(Na * Nb * Nc, 0.0), Py(Na * Nb * Nc, 0.0), Pz(Na * Nb * Nc, 0.0) {}

  // Spread a dipole vector at fractional grid position (wa, wb, wc)
  // onto the grid using B-spline weights.
  void spread(double mu_x, double mu_y, double mu_z,
              double wa, double wb, double wc) {
    const auto wts_a = bspline_weights(wa, order);
    const auto wts_b = bspline_weights(wb, order);
    const auto wts_c = bspline_weights(wc, order);
    const int i0a = static_cast<int>(std::floor(wa));
    const int i0b = static_cast<int>(std::floor(wb));
    const int i0c = static_cast<int>(std::floor(wc));
    const int iNa = static_cast<int>(Na);
    const int iNb = static_cast<int>(Nb);
    const int iNc = static_cast<int>(Nc);

    for (int ka = 0; ka < order; ++ka) {
      const int ia = ((i0a - ka) % iNa + iNa) % iNa;
      const double fa = wts_a[ka];
      for (int kb = 0; kb < order; ++kb) {
        const int ib = ((i0b - kb) % iNb + iNb) % iNb;
        const double fab = fa * wts_b[kb];
        for (int kc = 0; kc < order; ++kc) {
          const int ic = ((i0c - kc) % iNc + iNc) % iNc;
          const double w = fab * wts_c[kc];
          const size_t idx = ia * Nb * Nc + ib * Nc + ic;
          Px[idx] += mu_x * w;
          Py[idx] += mu_y * w;
          Pz[idx] += mu_z * w;
        }
      }
    }
  }

  void accumulate(const std::vector<trajan::core::EnhancedMolecule> &molecules,
                  const UnitCell &uc) {
    for (const auto &mol : molecules) {
      const occ::Vec3 com = mol.center_of_mass();

      // Compute dipole relative to COM (origin-independent for neutral molecules)
      double mu_x = 0.0, mu_y = 0.0, mu_z = 0.0;
      for (const auto &atom : mol.enhanced_atoms) {
        if (atom.charge == 0.0) continue;
        mu_x += atom.charge * (atom.x - com[0]);
        mu_y += atom.charge * (atom.y - com[1]);
        mu_z += atom.charge * (atom.z - com[2]);
      }

      // Fractional coordinate of COM, wrapped into [0, 1)
      occ::Vec3 frac = uc.to_fractional(com);
      frac = frac.array() - frac.array().floor();

      const double wa = frac[0] * static_cast<double>(Na);
      const double wb = frac[1] * static_cast<double>(Nb);
      const double wc = frac[2] * static_cast<double>(Nc);

      spread(mu_x, mu_y, mu_z, wa, wb, wc);
    }
    ++frame_count;
  }

  void write(const std::string &prefix, const UnitCell &uc) const {
    if (frame_count == 0) return;
    const double inv_n = 1.0 / static_cast<double>(frame_count);

    // Build averaged arrays
    const size_t N = Na * Nb * Nc;
    std::vector<double> avg_x(N), avg_y(N), avg_z(N), avg_mag(N);
    for (size_t i = 0; i < N; ++i) {
      avg_x[i] = Px[i] * inv_n;
      avg_y[i] = Py[i] * inv_n;
      avg_z[i] = Pz[i] * inv_n;
      avg_mag[i] = std::sqrt(avg_x[i]*avg_x[i] + avg_y[i]*avg_y[i] + avg_z[i]*avg_z[i]);
    }

    write_dipole_dx(prefix + "_Px.dx", avg_x, Na, Nb, Nc, uc,
                    "Time-averaged polarization Px (e*A)");
    write_dipole_dx(prefix + "_Py.dx", avg_y, Na, Nb, Nc, uc,
                    "Time-averaged polarization Py (e*A)");
    write_dipole_dx(prefix + "_Pz.dx", avg_z, Na, Nb, Nc, uc,
                    "Time-averaged polarization Pz (e*A)");
    write_dipole_dx(prefix + "_Pmag.dx", avg_mag, Na, Nb, Nc, uc,
                    "Time-averaged polarization magnitude |P| (e*A)");
  }
};

// ── Subcommand ─────────────────────────────────────────────────────────────

void run_dipole_subcommand(const DipoleOpts &opts, Trajectory &traj,
                           const Pipeline &pipeline) {
  // Validate grid options
  const bool do_grid = !opts.grid_out.empty();
  if (do_grid && opts.sel.empty())
    throw std::invalid_argument(
        "dipole: --sel is required when --grid-out is specified");
  if (do_grid && opts.order < 1)
    throw std::invalid_argument("dipole: --order must be >= 1");

  // Parse and validate molecule selector
  std::vector<trajan::io::SelectionCriteria> parsed_sel;
  if (!opts.sel.empty()) {
    parsed_sel = trajan::io::selection_validator(opts.sel,
                                                 trajan::core::ATOM_RESTRICTIONS);
  }

  size_t frame_count = 0;

  // Accumulators for whole-cell dipole
  double sum_dx{0}, sum_dy{0}, sum_dz{0};
  double sum_da{0}, sum_db{0}, sum_dc{0};
  double sum_total_charge{0};

  std::ofstream outfile;
  if (!opts.outfile.empty()) {
    outfile.open(opts.outfile);
    if (!outfile.is_open())
      throw std::runtime_error("dipole: cannot open output file: " + opts.outfile);
    outfile << "# frame  dx(e.A)  dy(e.A)  dz(e.A)  da(e.A)  db(e.A)  dc(e.A)\n";
  }

  // Grid setup — deferred until first frame (need unit cell)
  std::unique_ptr<DipoleGrid> grid;
  bool grid_initialised = false;
  bool charges_warned = false;

  trajan::log::Progress progress("dipole: reading frames");

  while (traj.next_frame()) {
    ++frame_count;
    progress.increment();
    pipeline.apply(traj.frame());

    const auto &atoms = traj.atoms();
    const auto &uc_opt = traj.unit_cell();

    // Warn once if charges are all zero
    if (!charges_warned) {
      bool any_nonzero = false;
      for (const auto &atom : atoms)
        if (atom.charge != 0.0) { any_nonzero = true; break; }
      if (!any_nonzero)
        trajan::log::warn("dipole: all atom charges are zero — use 'modify "
                          "--set-charge' to assign partial charges");
      charges_warned = true;
    }

    // ── Whole-cell dipole ────────────────────────────────────────────────
    double dx = 0.0, dy = 0.0, dz = 0.0, total_charge = 0.0;
    for (const auto &atom : atoms) {
      dx += atom.charge * atom.x;
      dy += atom.charge * atom.y;
      dz += atom.charge * atom.z;
      total_charge += atom.charge;
    }

    double da = dx, db = dy, dc = dz;
    if (uc_opt.has_value()) {
      const auto &direct = uc_opt->direct();
      const occ::Vec3 a_hat = direct.col(0).normalized();
      const occ::Vec3 b_hat = direct.col(1).normalized();
      const occ::Vec3 c_hat = direct.col(2).normalized();
      const occ::Vec3 d_vec(dx, dy, dz);
      da = d_vec.dot(a_hat);
      db = d_vec.dot(b_hat);
      dc = d_vec.dot(c_hat);
    } else if (frame_count == 1) {
      trajan::log::warn("dipole: no unit cell — a,b,c components equal x,y,z");
    }

    sum_dx += dx; sum_dy += dy; sum_dz += dz;
    sum_da += da; sum_db += db; sum_dc += dc;
    sum_total_charge += total_charge;

    if (outfile.is_open())
      outfile << fmt::format("{:6d}  {:12.6f}  {:12.6f}  {:12.6f}  {:12.6f}  "
                             "{:12.6f}  {:12.6f}\n",
                             frame_count - 1, dx, dy, dz, da, db, dc);

    // ── Polarization grid ────────────────────────────────────────────────
    if (do_grid) {
      if (!uc_opt.has_value())
        throw std::runtime_error("dipole --grid-out: unit cell required (CRYST1 record)");

      const UnitCell &uc = uc_opt.value();

      // Initialise grid on first frame
      if (!grid_initialised) {
        size_t Na, Nb, Nc;
        if (opts.grid.size() == 3) {
          Na = static_cast<size_t>(opts.grid[0]);
          Nb = static_cast<size_t>(opts.grid[1]);
          Nc = static_cast<size_t>(opts.grid[2]);
        } else {
          Na = fft_friendly(static_cast<size_t>(std::ceil(uc.a() / opts.spacing)));
          Nb = fft_friendly(static_cast<size_t>(std::ceil(uc.b() / opts.spacing)));
          Nc = fft_friendly(static_cast<size_t>(std::ceil(uc.c() / opts.spacing)));
        }
        trajan::log::info("dipole: polarization grid {}×{}×{} "
                          "(spacing ~{:.3f}/{:.3f}/{:.3f} Å), order={}",
                          Na, Nb, Nc,
                          uc.a() / Na, uc.b() / Nb, uc.c() / Nc,
                          opts.order);
        grid = std::make_unique<DipoleGrid>(Na, Nb, Nc, opts.order);
        grid_initialised = true;
      }

      const auto molecules = traj.get_molecules(parsed_sel);
      grid->accumulate(molecules, uc);
    }
  }

  progress.finish(fmt::format("dipole: processed {} frame(s)", frame_count));

  if (frame_count == 0) {
    trajan::log::warn("dipole: no frames available");
    return;
  }

  // ── Print whole-cell average ─────────────────────────────────────────────
  const double n = static_cast<double>(frame_count);
  const double avg_dx = sum_dx / n, avg_dy = sum_dy / n, avg_dz = sum_dz / n;
  const double avg_da = sum_da / n, avg_db = sum_db / n, avg_dc = sum_dc / n;
  const double mag = std::sqrt(avg_dx*avg_dx + avg_dy*avg_dy + avg_dz*avg_dz);

  constexpr size_t fw = 30;
  auto field = [&](const std::string &label, const std::string &value) {
    const int fill = static_cast<int>(fw) - static_cast<int>(label.size()) - 1;
    trajan::log::info("  |  {} {} : {}", label,
                      fill > 0 ? std::string(fill, '.') : "", value);
  };

  trajan::log::info("");
  trajan::log::info("  +- Dipole (average over {} frames) {}",
                    frame_count, std::string(38, '-'));
  field("total charge", fmt::format("{:+.6f} e", sum_total_charge / n));
  trajan::log::info("  |");
  field("dx", fmt::format("{:+.6f} e·Å", avg_dx));
  field("dy", fmt::format("{:+.6f} e·Å", avg_dy));
  field("dz", fmt::format("{:+.6f} e·Å", avg_dz));
  field("|d| (Cartesian)",
        fmt::format("{:.6f} e·Å  ({:.4f} D)", mag, mag * 4.80320));
  trajan::log::info("  |");
  field("da", fmt::format("{:+.6f} e·Å", avg_da));
  field("db", fmt::format("{:+.6f} e·Å", avg_db));
  field("dc", fmt::format("{:+.6f} e·Å", avg_dc));
  trajan::log::info("  +{}", std::string(77, '-'));

  if (outfile.is_open())
    trajan::log::info("dipole: per-frame data written to '{}'", opts.outfile);

  // ── Write polarization DX files ──────────────────────────────────────────
  if (do_grid && grid) {
    grid->write(opts.grid_out, traj.unit_cell().value());
    trajan::log::info("dipole: polarization field written to {}_P{{x,y,z,mag}}.dx",
                      opts.grid_out);
  }
}

CLI::App *add_dipole_subcommand(CLI::App &app, Trajectory &traj,
                                Pipeline &pipeline) {
  auto *dipole = app.add_subcommand(
      "dipole",
      "Compute the electric dipole moment of the simulation cell. "
      "With --grid-out, also outputs the time-averaged polarization vector "
      "field P(r) as OpenDX files. Requires load subcommand.");
  auto opts = std::make_shared<DipoleOpts>();

  dipole->add_option("-o,--output", opts->outfile,
                     "Write per-frame whole-cell dipole components to this file.");
  dipole->add_option("--sel", opts->sel,
                     "Molecule selector for polarization grid (molecule type or "
                     "index, e.g. 'mWAT'). Required with --grid-out.");
  dipole->add_option("--grid-out", opts->grid_out,
                     "Output prefix for polarization DX files. Produces "
                     "<prefix>_Px.dx, _Py.dx, _Pz.dx, _Pmag.dx.");
  dipole->add_option("--grid", opts->grid,
                     "Explicit grid dimensions Na Nb Nc (overrides --spacing).")
      ->expected(3);
  dipole->add_option("--spacing", opts->spacing,
                     "Target grid spacing in Å (derives Na Nb Nc from unit cell).")
      ->default_val(1.0)
      ->capture_default_str();
  dipole->add_option("--order", opts->order,
                     "B-spline assignment order for grid spreading (1=NGP, 4=cubic).")
      ->default_val(4)
      ->capture_default_str();

  dipole->callback([opts, &traj, &pipeline]() {
    trajan::log::set_subcommand_log_pattern("dipole");
    run_dipole_subcommand(*opts, traj, pipeline);
  });

  return dipole;
}

} // namespace trajan::main
