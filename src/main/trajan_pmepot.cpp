#include <trajan/core/frame.h>
#include <trajan/core/log.h>
#include <trajan/core/trajectory.h>
#include <trajan/main/trajan_pmepot.h>

#include <occ/3rdparty/pocketfft.h>
#include <occ/crystal/unitcell.h>

#include <cmath>
#include <complex>
#include <fstream>
#include <stdexcept>
#include <vector>

namespace trajan::main {

using occ::crystal::UnitCell;
using pocketfft::shape_t;
using pocketfft::stride_t;

static constexpr double PI = 3.14159265358979323846;
// Coulomb constant in practical units: V·Å / e
static constexpr double KE = 14.3996;

// ── B-spline helpers ──────────────────────────────────────────────────────

// Recursive B-spline M_p(x), support [0, p].
static double bspline(int p, double x) {
  if (x < 0.0 || x >= static_cast<double>(p))
    return 0.0;
  if (p == 1)
    return 1.0;
  return (x / (p - 1)) * bspline(p - 1, x) +
         ((p - x) / (p - 1)) * bspline(p - 1, x - 1.0);
}

// Compute p B-spline weights for a grid position w.
// Grid indices: (floor(w) - k) mod N, for k = 0..p-1.
// weight[k] = M_p(frac + k) where frac = w - floor(w).
static std::vector<double> bspline_weights(double w, int p) {
  int i0 = static_cast<int>(std::floor(w));
  double frac = w - i0;
  std::vector<double> wt(p);
  for (int k = 0; k < p; ++k)
    wt[k] = bspline(p, frac + k);
  return wt;
}

// |b(m, N, p)|^{-2}: reciprocal of the squared magnitude of the B-spline
// structure factor for grid axis with N points at reciprocal index m.
// b(m) = sum_{k=0}^{p-2} M_p(k+1) * exp(2πi·m·k/N)
static double bspline_inv_sq(int m, int N, int p) {
  const double theta = 2.0 * PI * m / N;
  std::complex<double> denom(0.0, 0.0);
  for (int k = 0; k < p - 1; ++k) {
    double wk = bspline(p, static_cast<double>(k + 1));
    denom += wk * std::polar(1.0, theta * k);
  }
  double mag2 = std::norm(denom);
  return (mag2 < 1e-20) ? 0.0 : 1.0 / mag2;
}

// ── Grid sizing ───────────────────────────────────────────────────────────

// Round n up to a 5-smooth number (only prime factors 2, 3, 5) so that
// pocketfft operates at peak efficiency.
static size_t fft_friendly(size_t n) {
  if (n < 1) n = 1;
  while (true) {
    size_t t = n;
    while (t % 2 == 0) t /= 2;
    while (t % 3 == 0) t /= 3;
    while (t % 5 == 0) t /= 5;
    if (t == 1)
      return n;
    ++n;
  }
}

// ── OpenDX output ─────────────────────────────────────────────────────────

static void write_dx(const std::string &path, const std::vector<double> &data,
                     size_t Na, size_t Nb, size_t Nc, const UnitCell &uc) {
  std::ofstream f(path);
  if (!f)
    throw std::runtime_error("pmepot: cannot open output file: " + path);

  // Grid step vectors in Cartesian Å (fractional step × lattice vector)
  const occ::Vec3 da = uc.a_vector() / static_cast<double>(Na);
  const occ::Vec3 db = uc.b_vector() / static_cast<double>(Nb);
  const occ::Vec3 dc = uc.c_vector() / static_cast<double>(Nc);

  f << "object 1 class gridpositions counts " << Na << " " << Nb << " " << Nc
    << "\n";
  f << "origin 0.0 0.0 0.0\n";
  f << "delta " << da[0] << " " << da[1] << " " << da[2] << "\n";
  f << "delta " << db[0] << " " << db[1] << " " << db[2] << "\n";
  f << "delta " << dc[0] << " " << dc[1] << " " << dc[2] << "\n";
  f << "object 2 class gridconnections counts " << Na << " " << Nb << " " << Nc
    << "\n";
  f << "object 3 class array type double rank 0 items " << (Na * Nb * Nc)
    << " data follows\n";

  size_t col = 0;
  for (const double v : data) {
    f << v;
    if (++col % 3 == 0)
      f << "\n";
    else
      f << "\t";
  }
  if (col % 3 != 0)
    f << "\n";

  f << "object \"electrostatic potential\" class field\n"
    << "component \"positions\" value 1\n"
    << "component \"connections\" value 2\n"
    << "component \"data\" value 3\n"
    << "attribute \"dep\" string \"positions\"\n";
}

// ── PmepotGrid ────────────────────────────────────────────────────────────

struct PmepotGrid {
  size_t Na, Nb, Nc, Nc2; // Nc2 = Nc/2+1 (r2c output size along c)
  int order;
  double beta;

  std::vector<float> Q;                    // charge mesh [Na·Nb·Nc]
  std::vector<std::complex<float>> Qk;     // reciprocal mesh [Na·Nb·Nc2]
  std::vector<float> green;                // Green's function [Na·Nb·Nc2]
  std::vector<double> accum;               // running potential sum [Na·Nb·Nc]
  size_t frame_count{0};

  PmepotGrid(size_t Na, size_t Nb, size_t Nc, int order, double beta)
      : Na(Na), Nb(Nb), Nc(Nc), Nc2(Nc / 2 + 1), order(order), beta(beta),
        Q(Na * Nb * Nc, 0.0f), Qk(Na * Nb * (Nc / 2 + 1)),
        green(Na * Nb * (Nc / 2 + 1), 0.0f), accum(Na * Nb * Nc, 0.0) {}

  // Build the reciprocal-space Green's function G(m1,m2,m3) once per cell.
  // G(m) = (4π·ke/V) · exp(-|k_m|²/4β²) / |k_m|² · |B(m1)|⁻²·|B(m2)|⁻²·|B(m3)|⁻²
  void build_green(const UnitCell &uc) {
    const double V = uc.volume();
    const double four_pi_ke_over_V = 4.0 * PI * KE / V;
    const double inv_4b2 = 1.0 / (4.0 * beta * beta);
    // Reciprocal matrix: column j is b*_j in Å⁻¹
    const occ::Mat3 recip = uc.reciprocal();

    // Precompute per-axis B-spline correction factors.
    std::vector<double> Ba(Na), Bb(Nb), Bc(Nc2);
    for (size_t m = 0; m < Na; ++m)
      Ba[m] = bspline_inv_sq(static_cast<int>(m), static_cast<int>(Na), order);
    for (size_t m = 0; m < Nb; ++m)
      Bb[m] = bspline_inv_sq(static_cast<int>(m), static_cast<int>(Nb), order);
    for (size_t m = 0; m < Nc2; ++m)
      Bc[m] = bspline_inv_sq(static_cast<int>(m), static_cast<int>(Nc), order);

    for (size_t ia = 0; ia < Na; ++ia) {
      // Map FFT index to signed reciprocal integer.
      const int m1 = (ia <= Na / 2) ? static_cast<int>(ia)
                                     : static_cast<int>(ia) - static_cast<int>(Na);
      for (size_t ib = 0; ib < Nb; ++ib) {
        const int m2 = (ib <= Nb / 2) ? static_cast<int>(ib)
                                       : static_cast<int>(ib) - static_cast<int>(Nb);
        for (size_t ic = 0; ic < Nc2; ++ic) {
          const int m3 = static_cast<int>(ic); // r2c: only non-negative m3

          const size_t idx = ia * Nb * Nc2 + ib * Nc2 + ic;

          if (m1 == 0 && m2 == 0 && m3 == 0) {
            green[idx] = 0.0f; // exclude DC (charge neutrality)
            continue;
          }

          const occ::Vec3 km =
              recip.col(0) * m1 + recip.col(1) * m2 + recip.col(2) * m3;
          const double k2 = km.squaredNorm();

          const double g = four_pi_ke_over_V * std::exp(-k2 * inv_4b2) / k2 *
                           Ba[ia] * Bb[ib] * Bc[ic];
          green[idx] = static_cast<float>(g);
        }
      }
    }
  }

  // Spread atomic partial charges onto Q[] using B-spline interpolation.
  void spread_charges(const std::vector<trajan::core::EnhancedAtom> &atoms,
                      const UnitCell &uc) {
    std::fill(Q.begin(), Q.end(), 0.0f);
    const int iNa = static_cast<int>(Na);
    const int iNb = static_cast<int>(Nb);
    const int iNc = static_cast<int>(Nc);

    for (const auto &atom : atoms) {
      if (atom.charge == 0.0)
        continue;

      occ::Vec3 frac = uc.to_fractional(occ::Vec3(atom.x, atom.y, atom.z));
      // Wrap into [0,1) along each axis.
      frac = frac.array() - frac.array().floor();

      const double wa = frac[0] * Na;
      const double wb = frac[1] * Nb;
      const double wc = frac[2] * Nc;

      const auto wts_a = bspline_weights(wa, order);
      const auto wts_b = bspline_weights(wb, order);
      const auto wts_c = bspline_weights(wc, order);

      const int i0a = static_cast<int>(std::floor(wa));
      const int i0b = static_cast<int>(std::floor(wb));
      const int i0c = static_cast<int>(std::floor(wc));

      for (int ka = 0; ka < order; ++ka) {
        const int ia = ((i0a - ka) % iNa + iNa) % iNa;
        const float wa_q = static_cast<float>(atom.charge * wts_a[ka]);
        for (int kb = 0; kb < order; ++kb) {
          const int ib = ((i0b - kb) % iNb + iNb) % iNb;
          const float wab_q = wa_q * static_cast<float>(wts_b[kb]);
          for (int kc = 0; kc < order; ++kc) {
            const int ic = ((i0c - kc) % iNc + iNc) % iNc;
            Q[ia * iNb * iNc + ib * iNc + ic] +=
                wab_q * static_cast<float>(wts_c[kc]);
          }
        }
      }
    }
  }

  // Forward r2c FFT → apply Green's function → inverse c2r → accumulate.
  void fft_and_accumulate() {
    const std::ptrdiff_t sf = sizeof(float);
    const std::ptrdiff_t sc = sizeof(std::complex<float>);

    const shape_t shape_r = {Na, Nb, Nc};
    const stride_t stride_r = {static_cast<std::ptrdiff_t>(Nb * Nc) * sf,
                                static_cast<std::ptrdiff_t>(Nc) * sf, sf};
    const stride_t stride_c = {static_cast<std::ptrdiff_t>(Nb * Nc2) * sc,
                                static_cast<std::ptrdiff_t>(Nc2) * sc, sc};

    // Forward real-to-complex along all three axes (last axis r2c, rest c2c).
    pocketfft::r2c(shape_r, stride_r, stride_c, {0, 1, 2}, true,
                   Q.data(), Qk.data(), 1.0f);

    // Multiply by Green's function in reciprocal space.
    for (size_t i = 0; i < Qk.size(); ++i)
      Qk[i] *= green[i];

    // Inverse complex-to-real, with 1/N normalisation.
    const float fct = 1.0f / static_cast<float>(Na * Nb * Nc);
    std::vector<float> phi(Na * Nb * Nc);
    pocketfft::c2r(shape_r, stride_c, stride_r, {0, 1, 2}, false,
                   Qk.data(), phi.data(), fct);

    // Accumulate in double precision.
    for (size_t i = 0; i < phi.size(); ++i)
      accum[i] += static_cast<double>(phi[i]);

    ++frame_count;
  }

  std::vector<double> average() const {
    std::vector<double> avg(accum.size());
    const double inv_n = frame_count > 0 ? 1.0 / frame_count : 1.0;
    for (size_t i = 0; i < accum.size(); ++i)
      avg[i] = accum[i] * inv_n;
    return avg;
  }
};

// ── Subcommand ────────────────────────────────────────────────────────────

void run_pmepot_subcommand(const PmepotOpts &opts, Trajectory &traj,
                           const Pipeline &pipeline) {
  if (!traj.next_frame())
    throw std::runtime_error("pmepot: no frames available");

  trajan::core::Frame &frame = traj.frame();
  if (!frame.has_unit_cell())
    throw std::runtime_error(
        "pmepot: frame has no unit cell (CRYST1 record required)");

  const UnitCell &uc0 = frame.unit_cell().value();

  // Determine grid dimensions — fixed for the whole trajectory.
  size_t Na, Nb, Nc;
  if (opts.grid.size() == 3) {
    Na = static_cast<size_t>(opts.grid[0]);
    Nb = static_cast<size_t>(opts.grid[1]);
    Nc = static_cast<size_t>(opts.grid[2]);
  } else {
    Na = fft_friendly(static_cast<size_t>(std::ceil(uc0.a() / opts.spacing)));
    Nb = fft_friendly(static_cast<size_t>(std::ceil(uc0.b() / opts.spacing)));
    Nc = fft_friendly(static_cast<size_t>(std::ceil(uc0.c() / opts.spacing)));
  }

  trajan::log::info(
      "pmepot: grid {}×{}×{} (spacing ~{:.3f}/{:.3f}/{:.3f} Å), "
      "β={:.4f} Å⁻¹, order={}",
      Na, Nb, Nc, uc0.a() / Na, uc0.b() / Nb, uc0.c() / Nc,
      opts.ewaldcof, opts.order);

  PmepotGrid grid(Na, Nb, Nc, opts.order, opts.ewaldcof);
  grid.build_green(uc0);

  // Track the last cell to detect changes (NPT trajectories).
  UnitCell last_uc = uc0;

  auto process_frame = [&](trajan::core::Frame &fr) {
    const UnitCell &uc = fr.unit_cell().value();

    // Recompute Green's function only when the cell has changed noticeably.
    if (std::abs(uc.a() - last_uc.a()) > 1e-4 ||
        std::abs(uc.b() - last_uc.b()) > 1e-4 ||
        std::abs(uc.c() - last_uc.c()) > 1e-4 ||
        std::abs(uc.alpha() - last_uc.alpha()) > 1e-6 ||
        std::abs(uc.beta()  - last_uc.beta())  > 1e-6 ||
        std::abs(uc.gamma() - last_uc.gamma()) > 1e-6) {
      trajan::log::debug("pmepot: cell changed at frame {}, rebuilding Green's function",
                         grid.frame_count + 1);
      grid.build_green(uc);
      last_uc = uc;
    }

    grid.spread_charges(fr.atoms(), uc);
    grid.fft_and_accumulate();
  };

  pipeline.apply(frame);
  process_frame(frame);
  while (traj.next_frame()) {
    pipeline.apply(traj.frame());
    process_frame(traj.frame());
  }

  trajan::log::info("pmepot: averaged over {} frame(s)", grid.frame_count);

  const auto avg = grid.average();
  write_dx(opts.outfile, avg, Na, Nb, Nc, last_uc);
  trajan::log::info("pmepot: wrote potential to {}", opts.outfile);
}

CLI::App *add_pmepot_subcommand(CLI::App &app, Trajectory &traj,
                                Pipeline &pipeline) {
  CLI::App *cmd = app.add_subcommand(
      "pmepot",
      "Compute the time-averaged electrostatic potential on a 3-D grid via "
      "the reciprocal Ewald sum (PMEPOT). Partial charges must be set on "
      "atoms (modify --set-charge). Output is an OpenDX file (.dx) in Volts. "
      "Requires load subcommand.");

  auto opts = std::make_shared<PmepotOpts>();
  cmd->add_option("outfile", opts->outfile, "Output OpenDX file (.dx)")
      ->required();
  cmd->add_option("--spacing", opts->spacing,
                  "Target grid spacing in Å (derives Na Nb Nc from unit cell)")
      ->default_val(1.0)
      ->capture_default_str();
  cmd->add_option("--grid", opts->grid,
                  "Explicit grid dimensions Na Nb Nc along a, b, c "
                  "(overrides --spacing)")
      ->expected(3);
  cmd->add_option("--ewaldcof,--beta", opts->ewaldcof,
                  "Ewald coefficient β in Å⁻¹ (match the value used in the "
                  "MD simulation; typical PME value: 0.25–0.35)")
      ->default_val(0.25)
      ->capture_default_str();
  cmd->add_option("--order", opts->order,
                  "B-spline interpolation order (4 = cubic, must be ≥ 2)")
      ->default_val(4)
      ->capture_default_str();

  cmd->callback([opts, &traj, &pipeline]() {
    trajan::log::set_subcommand_log_pattern("pmepot");
    trajan::log::debug("Beginning pmepot subcommand");
    if (opts->order < 2)
      throw std::invalid_argument("pmepot: --order must be >= 2");
    run_pmepot_subcommand(*opts, traj, pipeline);
  });

  return cmd;
}

} // namespace trajan::main
