#include <trajan/core/frame.h>
#include <trajan/core/log.h>
#include <trajan/core/trajectory.h>
#include <trajan/main/trajan_pmepot.h>

#include <occ/3rdparty/pocketfft.h>
#include <occ/core/constants.h>
#include <occ/crystal/unitcell.h>

#include <tbb/combinable.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <fstream>
#include <sstream>
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

  f << "# PME electrostatic potential (Volts)\n";
  f << "object 1 class gridpositions counts " << Na << " " << Nb << " " << Nc
    << "\n";
  // VMD centers the grid at Cartesian (0,0,0): origin = -((Na-1)/2*da + (Nb-1)/2*db + (Nc-1)/2*dc)
  const occ::Vec3 origin = -(da * ((Na - 1) * 0.5) + db * ((Nb - 1) * 0.5) + dc * ((Nc - 1) * 0.5));
  f << "origin " << origin[0] << " " << origin[1] << " " << origin[2] << "\n";
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
    // OCC reciprocal() uses the crystallographic convention (no 2π factor).
    // PME requires physics wavevectors k = 2π·g, so multiply by 2π.
    const occ::Mat3 recip = uc.reciprocal() * occ::constants::two_pi<double>;

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

      // Shift by (N-1)/2 so that grid index 0 corresponds to the VMD origin.
      // This matches VMD's convention: origin = -((Na-1)/2*da + ...).
      const double wa = frac[0] * Na + (Na - 1) * 0.5;
      const double wb = frac[1] * Nb + (Nb - 1) * 0.5;
      const double wc = frac[2] * Nc + (Nc - 1) * 0.5;

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

    // G already carries 1/V; the IDFT should not apply an additional 1/N.
    // Using fct=1 gives the unnormalized sum that matches the continuous formula
    // φ(r) = (1/V) Σ_m G(k_m)·Q̃(m)·exp(ik_m·r) where G = 4πke/V·exp(…)/k².
    const float fct = 1.0f;
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

// ── RealspacePotGrid ──────────────────────────────────────────────────────

struct RealspacePotGrid {
  size_t Na, Nb, Nc;
  double cutoff;
  std::vector<double> accum;
  size_t frame_count{0};

  RealspacePotGrid(size_t Na, size_t Nb, size_t Nc, double cutoff)
      : Na(Na), Nb(Nb), Nc(Nc), cutoff(cutoff), accum(Na * Nb * Nc, 0.0) {}

  // For each charged atom, directly search the nearby grid-index box and
  // accumulate KE·q/r. Parallelised over atoms via TBB; thread-local
  // accumulation avoids races without atomics.
  void accumulate_frame(const std::vector<trajan::core::EnhancedAtom> &atoms,
                        const UnitCell &uc) {
    std::vector<double> charges;
    std::vector<occ::Vec3> atom_carts;
    for (const auto &a : atoms) {
      if (a.charge == 0.0) continue;
      charges.push_back(a.charge);
      atom_carts.emplace_back(a.x, a.y, a.z);
    }
    if (charges.empty()) { ++frame_count; return; }

    const size_t N_charged = charges.size();
    trajan::log::debug("pmepot: frame {}: {} charged atoms, {}×{}×{}={} grid points",
                       frame_count + 1, N_charged, Na, Nb, Nc, Na * Nb * Nc);

    const occ::Mat3 D    = uc.direct();
    const occ::Mat3 Dinv = D.inverse();
    // Search radius in grid cells along each lattice direction (conservative).
    const int ra  = static_cast<int>(std::ceil(cutoff * Na / uc.a())) + 1;
    const int rb  = static_cast<int>(std::ceil(cutoff * Nb / uc.b())) + 1;
    const int rc  = static_cast<int>(std::ceil(cutoff * Nc / uc.c())) + 1;
    const int iNa = static_cast<int>(Na);
    const int iNb = static_cast<int>(Nb);
    const int iNc = static_cast<int>(Nc);
    const double off_a   = (Na - 1) * 0.5;
    const double off_b   = (Nb - 1) * 0.5;
    const double off_c   = (Nc - 1) * 0.5;
    const double cutoffsq = cutoff * cutoff;
    trajan::log::debug("pmepot: search box ±{}×{}×{} cells, starting TBB loop", ra, rb, rc);

    auto t0 = std::chrono::steady_clock::now();

    // Each TBB thread accumulates into its own private vector to avoid races.
    tbb::combinable<std::vector<double>> local_accum(
        [&] { return std::vector<double>(Na * Nb * Nc, 0.0); });

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, N_charged),
        [&](const tbb::blocked_range<size_t> &range) {
          auto &lacc = local_accum.local();
          for (size_t j = range.begin(); j < range.end(); ++j) {
            const occ::Vec3 &ac = atom_carts[j];
            // Atom fractional position, wrapped to [0,1).
            occ::Vec3 af = Dinv * ac;
            af.array() -= af.array().floor();
            // Grid index nearest to the atom's fractional position.
            const int ci_a = static_cast<int>(std::round(af[0] * Na + off_a));
            const int ci_b = static_cast<int>(std::round(af[1] * Nb + off_b));
            const int ci_c = static_cast<int>(std::round(af[2] * Nc + off_c));

            for (int da = -ra; da <= ra; ++da) {
              const int ia = ((ci_a + da) % iNa + iNa) % iNa;
              // Fractional difference along a (minimum image applied below).
              const double dfa = (ia - off_a) / Na - af[0];
              for (int db = -rb; db <= rb; ++db) {
                const int ib = ((ci_b + db) % iNb + iNb) % iNb;
                const double dfb = (ib - off_b) / Nb - af[1];
                for (int dc = -rc; dc <= rc; ++dc) {
                  const int ic = ((ci_c + dc) % iNc + iNc) % iNc;
                  const double dfc = (ic - off_c) / Nc - af[2];
                  // Minimum image in fractional space, then to Cartesian.
                  const occ::Vec3 df(dfa - std::round(dfa),
                                     dfb - std::round(dfb),
                                     dfc - std::round(dfc));
                  const double rsq = (D * df).squaredNorm();
                  if (rsq < 1e-10 || rsq > cutoffsq) continue;
                  lacc[ia * Nb * Nc + ib * Nc + ic] +=
                      KE * charges[j] / std::sqrt(rsq);
                }
              }
            }
          }
        });

    // Reduce thread-local results into the running accumulator.
    local_accum.combine_each([&](const std::vector<double> &la) {
      for (size_t i = 0; i < Na * Nb * Nc; ++i)
        accum[i] += la[i];
    });

    auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(
                  std::chrono::steady_clock::now() - t0).count();
    trajan::log::debug("pmepot: accumulate_frame done in {}ms", dt);
    ++frame_count;
  }

  // pmepot_custom style: spread charges onto Q via B-splines, then for each
  // grid point sum KE·Q[neighbor]/dist over all neighbors within ±cell_radius
  // cells. The inv_dist table is precomputed once per frame (one sqrt per
  // unique offset), so the inner loop is only a multiply-add per neighbor.
  void accumulate_frame_grid_spread(
      const std::vector<trajan::core::EnhancedAtom> &atoms,
      const UnitCell &uc, int cell_radius, int order) {

    trajan::log::debug("pmepot: frame {}: grid-spread mode, cell_radius={}",
                       frame_count + 1, cell_radius);
    auto t0 = std::chrono::steady_clock::now();

    // ── 1. Spread atomic charges onto Q via B-splines (VMD convention) ──
    std::vector<float> Q(Na * Nb * Nc, 0.0f);
    const int iNa = static_cast<int>(Na);
    const int iNb = static_cast<int>(Nb);
    const int iNc = static_cast<int>(Nc);
    for (const auto &atom : atoms) {
      if (atom.charge == 0.0) continue;
      occ::Vec3 frac = uc.to_fractional(occ::Vec3(atom.x, atom.y, atom.z));
      frac = frac.array() - frac.array().floor();
      const double wa = frac[0] * Na + (Na - 1) * 0.5;
      const double wb = frac[1] * Nb + (Nb - 1) * 0.5;
      const double wc = frac[2] * Nc + (Nc - 1) * 0.5;
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
    auto t1 = std::chrono::steady_clock::now();
    trajan::log::debug("pmepot: charge spreading done in {}ms",
                       std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count());

    // ── 2. Precompute inv_dist for every (da,db,dc) offset in [-R,R]³ ──
    // The Cartesian vector from grid point (ia,ib,ic) to its neighbor at
    // offset (da,db,dc) is D·(da/Na, db/Nb, dc/Nc) — constant for fixed offset.
    const int R   = cell_radius;
    const int box = 2 * R + 1;
    const occ::Mat3 D = uc.direct();
    std::vector<double> inv_dist(box * box * box, 0.0);
    for (int da = -R; da <= R; ++da)
      for (int db = -R; db <= R; ++db)
        for (int dc = -R; dc <= R; ++dc) {
          if (da == 0 && db == 0 && dc == 0) continue;
          const occ::Vec3 dr = D * occ::Vec3(static_cast<double>(da) / Na,
                                             static_cast<double>(db) / Nb,
                                             static_cast<double>(dc) / Nc);
          inv_dist[(da + R) * box * box + (db + R) * box + (dc + R)] =
              1.0 / dr.norm();
        }

    // ── 3. Convolve Q with inv_dist, accumulate into accum ──
    // Parallelise over the outer (ia) dimension — no race since each output
    // grid point maps to a unique accum index.
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, Na),
        [&](const tbb::blocked_range<size_t> &range) {
          for (size_t ia = range.begin(); ia < range.end(); ++ia) {
            for (size_t ib = 0; ib < Nb; ++ib) {
              for (size_t ic = 0; ic < Nc; ++ic) {
                double phi = 0.0;
                for (int da = -R; da <= R; ++da) {
                  const int ia2 = (static_cast<int>(ia) + da + iNa * (R + 1)) % iNa;
                  for (int db = -R; db <= R; ++db) {
                    const int ib2 = (static_cast<int>(ib) + db + iNb * (R + 1)) % iNb;
                    for (int dc = -R; dc <= R; ++dc) {
                      if (da == 0 && db == 0 && dc == 0) continue;
                      const int ic2 = (static_cast<int>(ic) + dc + iNc * (R + 1)) % iNc;
                      phi += Q[ia2 * iNb * iNc + ib2 * iNc + ic2] *
                             inv_dist[(da + R) * box * box + (db + R) * box + (dc + R)];
                    }
                  }
                }
                accum[ia * Nb * Nc + ib * Nc + ic] += KE * phi;
              }
            }
          }
        });

    auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(
                  std::chrono::steady_clock::now() - t1).count();
    trajan::log::debug("pmepot: grid convolution done in {}ms", dt);
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

  if (opts.real_space && opts.grid_spread) {
    trajan::log::info(
        "pmepot: grid-spread real-space on {}×{}×{} grid "
        "(spacing ~{:.3f}/{:.3f}/{:.3f} Å), B-spline order={}, cell_radius={} ({}×{}×{} neighborhood)",
        Na, Nb, Nc, uc0.a() / Na, uc0.b() / Nb, uc0.c() / Nc,
        opts.order, opts.cell_radius,
        2 * opts.cell_radius + 1, 2 * opts.cell_radius + 1, 2 * opts.cell_radius + 1);
  } else if (opts.real_space) {
    trajan::log::info(
        "pmepot: real-space direct Coulomb sum on {}×{}×{} grid "
        "(spacing ~{:.3f}/{:.3f}/{:.3f} Å), cutoff={:.2f} Å",
        Na, Nb, Nc, uc0.a() / Na, uc0.b() / Nb, uc0.c() / Nc, opts.cutoff);
  } else {
    trajan::log::info(
        "pmepot: reciprocal-space PME on {}×{}×{} grid "
        "(spacing ~{:.3f}/{:.3f}/{:.3f} Å), β={:.4f} Å⁻¹, order={}",
        Na, Nb, Nc, uc0.a() / Na, uc0.b() / Nb, uc0.c() / Nc,
        opts.ewaldcof, opts.order);
  }

  UnitCell last_uc = uc0;
  trajan::log::Progress progress("pmepot: frames processed");

  if (opts.real_space) {
    RealspacePotGrid grid(Na, Nb, Nc, opts.cutoff);

    auto process_frame = [&](trajan::core::Frame &fr) {
      if (opts.grid_spread)
        grid.accumulate_frame_grid_spread(fr.atoms(), fr.unit_cell().value(),
                                          opts.cell_radius, opts.order);
      else
        grid.accumulate_frame(fr.atoms(), fr.unit_cell().value());
    };

    pipeline.apply(frame, "pmepot");
    process_frame(frame);
    if (grid.frame_count % 100 == 0)
      progress.update(static_cast<int64_t>(grid.frame_count));

    while (traj.next_frame()) {
      pipeline.apply(traj.frame(), "pmepot");
      process_frame(traj.frame());
      last_uc = traj.frame().unit_cell().value();
      if (grid.frame_count % 100 == 0)
        progress.update(static_cast<int64_t>(grid.frame_count));
    }

    progress.finish(fmt::format("pmepot: averaged over {} frame(s)", grid.frame_count));
    trajan::log::info("pmepot: averaged over {} frame(s)", grid.frame_count);
    write_dx(opts.outfile, grid.average(), Na, Nb, Nc, last_uc);
  } else {
    PmepotGrid grid(Na, Nb, Nc, opts.order, opts.ewaldcof);
    grid.build_green(uc0);

    auto process_frame = [&](trajan::core::Frame &fr) {
      const UnitCell &uc = fr.unit_cell().value();
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

    pipeline.apply(frame, "pmepot");
    process_frame(frame);
    if (grid.frame_count % 100 == 0)
      progress.update(static_cast<int64_t>(grid.frame_count));

    while (traj.next_frame()) {
      pipeline.apply(traj.frame(), "pmepot");
      process_frame(traj.frame());
      if (grid.frame_count % 100 == 0)
        progress.update(static_cast<int64_t>(grid.frame_count));
    }

    progress.finish(fmt::format("pmepot: averaged over {} frame(s)", grid.frame_count));
    trajan::log::info("pmepot: averaged over {} frame(s)", grid.frame_count);
    write_dx(opts.outfile, grid.average(), Na, Nb, Nc, last_uc);
  }

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
  cmd->add_flag("--real-space", opts->real_space,
                "Compute potential by direct real-space Coulomb sum (O(N·Ng)) "
                "instead of PME.");
  cmd->add_option("--cutoff", opts->cutoff,
                  "Neighbour cutoff radius in Å for --real-space mode.")
      ->default_val(10.0)
      ->capture_default_str();
  cmd->add_flag("--grid-spread", opts->grid_spread,
                "Grid-space real-space (pmepot_custom style): spread charges "
                "onto the grid via B-splines, then sum KE·Q[neighbor]/dist "
                "over a fixed cell neighborhood. Requires --real-space. "
                "--cutoff is ignored; use --cell-radius instead.");
  cmd->add_option("--cell-radius", opts->cell_radius,
                  "Neighbor search radius in grid cells for --grid-spread "
                  "(default 2 = 5×5×5 = 125 neighbors, matching pmepot_custom).")
      ->default_val(2)
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

// ── DX reader ─────────────────────────────────────────────────────────────

struct DXGrid {
  size_t Na{0}, Nb{0}, Nc{0};
  occ::Vec3 origin{0, 0, 0};
  occ::Vec3 da{0, 0, 0}, db{0, 0, 0}, dc{0, 0, 0};
  std::vector<double> data; // Na*Nb*Nc row-major [ia][ib][ic]

  double at(size_t ia, size_t ib, size_t ic) const {
    return data[ia * Nb * Nc + ib * Nc + ic];
  }
};

static DXGrid read_dx(const std::string &path) {
  std::ifstream f(path);
  if (!f)
    throw std::runtime_error("dxreduce: cannot open file: " + path);

  DXGrid g;
  std::string line;
  int delta_count = 0;
  bool reading_data = false;
  size_t data_expected = 0, data_read = 0;

  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#')
      continue;

    if (!reading_data) {
      if (line.find("object 1 class gridpositions counts") != std::string::npos) {
        std::sscanf(line.c_str(),
                    "object 1 class gridpositions counts %zu %zu %zu",
                    &g.Na, &g.Nb, &g.Nc);
      } else if (line.substr(0, 6) == "origin") {
        std::sscanf(line.c_str(), "origin %lf %lf %lf",
                    &g.origin[0], &g.origin[1], &g.origin[2]);
      } else if (line.substr(0, 5) == "delta") {
        if (delta_count == 0)
          std::sscanf(line.c_str(), "delta %lf %lf %lf",
                      &g.da[0], &g.da[1], &g.da[2]);
        else if (delta_count == 1)
          std::sscanf(line.c_str(), "delta %lf %lf %lf",
                      &g.db[0], &g.db[1], &g.db[2]);
        else
          std::sscanf(line.c_str(), "delta %lf %lf %lf",
                      &g.dc[0], &g.dc[1], &g.dc[2]);
        ++delta_count;
      } else if (line.find("data follows") != std::string::npos) {
        auto p = line.find("items");
        if (p != std::string::npos)
          data_expected = std::stoul(line.substr(p + 6));
        if (data_expected == 0)
          data_expected = g.Na * g.Nb * g.Nc;
        g.data.resize(data_expected, 0.0);
        reading_data = true;
      }
    } else {
      std::istringstream iss(line);
      double v;
      while (iss >> v && data_read < data_expected)
        g.data[data_read++] = v;
      if (data_read >= data_expected)
        break;
    }
  }

  if (g.Na == 0 || g.Nb == 0 || g.Nc == 0)
    throw std::runtime_error("dxreduce: failed to read grid dimensions from " + path);
  if (data_read != data_expected)
    throw std::runtime_error(
        fmt::format("dxreduce: expected {} data values, got {}", data_expected, data_read));

  return g;
}

// ── Axis averaging & Cartesian rebinning ──────────────────────────────────

enum class AxisType { Lattice, Cartesian };

struct ParsedAxes {
  AxisType type;
  std::vector<int> indices; // 0=a/x, 1=b/y, 2=c/z
};

// Parse "a","bc","xyz" etc. Disallows mixing lattice and Cartesian labels.
static ParsedAxes parse_axes(const std::string &s) {
  bool has_cart = s.find_first_of("xyz") != std::string::npos;
  bool has_latt = s.find_first_of("abc") != std::string::npos;
  if (has_cart && has_latt)
    throw std::invalid_argument(
        "dxreduce: cannot mix a/b/c and x/y/z in --average");
  std::vector<int> axes;
  for (char c : s) {
    if (c == 'a' || c == 'x')      axes.push_back(0);
    else if (c == 'b' || c == 'y') axes.push_back(1);
    else if (c == 'c' || c == 'z') axes.push_back(2);
    else
      throw std::invalid_argument(
          fmt::format("dxreduce: unknown axis '{}' (use a,b,c or x,y,z)", c));
  }
  std::sort(axes.begin(), axes.end());
  axes.erase(std::unique(axes.begin(), axes.end()), axes.end());
  return {has_cart ? AxisType::Cartesian : AxisType::Lattice, axes};
}

// True when all three delta vectors are axis-aligned (orthogonal cell).
// In that case Cartesian x/y/z maps 1-to-1 onto lattice a/b/c.
static bool is_orthogonal(const DXGrid &g, double tol = 1e-6) {
  return std::abs(g.da[1]) < tol && std::abs(g.da[2]) < tol &&
         std::abs(g.db[0]) < tol && std::abs(g.db[2]) < tol &&
         std::abs(g.dc[0]) < tol && std::abs(g.dc[1]) < tol;
}

// Fast lattice-axis averaging (exact, no rebinning).
static std::pair<std::vector<double>, std::vector<int>>
average_lattice_axes(const DXGrid &g, const std::vector<int> &avg_axes) {
  std::vector<int> remain;
  for (int k = 0; k < 3; ++k)
    if (std::find(avg_axes.begin(), avg_axes.end(), k) == avg_axes.end())
      remain.push_back(k);

  const size_t Ns[3] = {g.Na, g.Nb, g.Nc};
  size_t out_size = 1;
  for (int r : remain) out_size *= Ns[r];

  double denom = 1.0;
  for (int k : avg_axes) denom *= static_cast<double>(Ns[k]);

  std::vector<double> result(out_size, 0.0);
  for (size_t ia = 0; ia < g.Na; ++ia) {
    for (size_t ib = 0; ib < g.Nb; ++ib) {
      for (size_t ic = 0; ic < g.Nc; ++ic) {
        const size_t idx[3] = {ia, ib, ic};
        size_t out_idx = 0;
        for (size_t r = 0; r < remain.size(); ++r) {
          size_t stride = 1;
          for (size_t s = r + 1; s < remain.size(); ++s)
            stride *= Ns[remain[s]];
          out_idx += idx[remain[r]] * stride;
        }
        result[out_idx] += g.at(ia, ib, ic);
      }
    }
  }
  for (double &v : result) v /= denom;
  return {result, remain};
}

// Cartesian rebinning for non-orthogonal cells.
// Projects each grid point onto the remaining Cartesian axis/axes, bins by
// the finest grid step along that direction, then averages within bins.
// Returns (values, pos0_vec, pos1_vec, N0, N1).
// For 1D: N1=1, pos1_vec is empty.
struct CartRebin {
  std::vector<double> values;
  std::vector<double> pos0, pos1;
  size_t N0{0}, N1{1};
};

static CartRebin rebin_cartesian(const DXGrid &g, const std::vector<int> &rem_cart) {
  const occ::Vec3 *deltas[3] = {&g.da, &g.db, &g.dc};
  const bool is_2d = (rem_cart.size() == 2);
  const int ax0 = rem_cart[0];
  const int ax1 = is_2d ? rem_cart[1] : -1;

  // Finest non-zero step contribution along a Cartesian axis.
  auto finest_step = [&](int cart_ax) {
    double step = std::numeric_limits<double>::max();
    for (int k = 0; k < 3; ++k) {
      double c = std::abs((*deltas[k])[cart_ax]);
      if (c > 1e-10) step = std::min(step, c);
    }
    return step;
  };

  auto cart_coord = [&](size_t ia, size_t ib, size_t ic, int cart_ax) {
    return ia * (*deltas[0])[cart_ax] +
           ib * (*deltas[1])[cart_ax] +
           ic * (*deltas[2])[cart_ax];
  };

  const double step0 = finest_step(ax0);
  const double step1 = is_2d ? finest_step(ax1) : 1.0;

  // Find coordinate ranges by visiting all grid points.
  double mn0 = std::numeric_limits<double>::max(), mx0 = std::numeric_limits<double>::lowest();
  double mn1 = std::numeric_limits<double>::max(), mx1 = std::numeric_limits<double>::lowest();
  for (size_t ia = 0; ia < g.Na; ++ia) {
    for (size_t ib = 0; ib < g.Nb; ++ib) {
      for (size_t ic = 0; ic < g.Nc; ++ic) {
        double p0 = cart_coord(ia, ib, ic, ax0);
        mn0 = std::min(mn0, p0); mx0 = std::max(mx0, p0);
        if (is_2d) {
          double p1 = cart_coord(ia, ib, ic, ax1);
          mn1 = std::min(mn1, p1); mx1 = std::max(mx1, p1);
        }
      }
    }
  }

  const size_t N0 = static_cast<size_t>(std::round((mx0 - mn0) / step0)) + 1;
  const size_t N1 = is_2d ? static_cast<size_t>(std::round((mx1 - mn1) / step1)) + 1 : 1;

  std::vector<double> sums(N0 * N1, 0.0);
  std::vector<size_t> counts(N0 * N1, 0);

  for (size_t ia = 0; ia < g.Na; ++ia) {
    for (size_t ib = 0; ib < g.Nb; ++ib) {
      for (size_t ic = 0; ic < g.Nc; ++ic) {
        double p0 = cart_coord(ia, ib, ic, ax0);
        size_t b0 = static_cast<size_t>(std::round((p0 - mn0) / step0));
        size_t b1 = 0;
        if (is_2d) {
          double p1 = cart_coord(ia, ib, ic, ax1);
          b1 = static_cast<size_t>(std::round((p1 - mn1) / step1));
        }
        if (b0 < N0 && b1 < N1) {
          sums[b0 * N1 + b1] += g.at(ia, ib, ic);
          ++counts[b0 * N1 + b1];
        }
      }
    }
  }

  std::vector<double> values(N0 * N1, 0.0);
  for (size_t i = 0; i < N0 * N1; ++i)
    values[i] = counts[i] > 0 ? sums[i] / static_cast<double>(counts[i]) : 0.0;

  std::vector<double> p0v(N0), p1v(is_2d ? N1 : 0);
  for (size_t i = 0; i < N0; ++i) p0v[i] = mn0 + i * step0;
  for (size_t i = 0; i < N1 && is_2d; ++i) p1v[i] = mn1 + i * step1;

  return {values, p0v, p1v, N0, N1};
}

// ── Output writers ─────────────────────────────────────────────────────────

// Generic 1D writer: explicit position vector, arbitrary axis label.
static void write_1d(const std::string &path,
                     const std::vector<double> &positions,
                     const std::vector<double> &values,
                     const std::string &pos_label,
                     const std::string &avg_label) {
  std::ofstream f(path);
  if (!f) throw std::runtime_error("dxreduce: cannot open output: " + path);
  f << fmt::format("# dxreduce 1D: remaining axis {}\n", pos_label);
  f << fmt::format("# averaged over: {}\n", avg_label);
  f << "# pos(A)  potential(V)\n";
  for (size_t i = 0; i < values.size(); ++i)
    f << fmt::format("{:12.6f}  {:16.8f}\n", positions[i], values[i]);
}

// Generic 2D writer: gnuplot pm3d compatible, blank line between i0 blocks.
static void write_2d(const std::string &path,
                     const std::vector<double> &pos0,
                     const std::vector<double> &pos1,
                     const std::vector<double> &values,
                     const std::string &label0,
                     const std::string &label1,
                     const std::string &avg_label) {
  std::ofstream f(path);
  if (!f) throw std::runtime_error("dxreduce: cannot open output: " + path);
  const size_t N0 = pos0.size(), N1 = pos1.size();
  f << fmt::format("# dxreduce 2D: axes {} {}\n", label0, label1);
  f << fmt::format("# averaged over: {}\n", avg_label);
  f << fmt::format("# N{}={}  N{}={}\n", label0, N0, label1, N1);
  f << fmt::format("# pos_{}(A)  pos_{}(A)  potential(V)\n", label0, label1);
  for (size_t i0 = 0; i0 < N0; ++i0) {
    for (size_t i1 = 0; i1 < N1; ++i1)
      f << fmt::format("{:12.6f}  {:12.6f}  {:16.8f}\n",
                       pos0[i0], pos1[i1], values[i0 * N1 + i1]);
    f << "\n";
  }
}

// ── Subcommand ─────────────────────────────────────────────────────────────

void run_dxreduce_subcommand(const DXReduceOpts &opts) {
  trajan::log::info("dxreduce: reading {}", opts.infile);
  const DXGrid g = read_dx(opts.infile);
  trajan::log::info("dxreduce: grid {}×{}×{}", g.Na, g.Nb, g.Nc);

  if (opts.average_axes.empty())
    throw std::invalid_argument("dxreduce: --average must specify at least one axis");

  const auto [atype, avg_idx] = parse_axes(opts.average_axes);
  if (avg_idx.size() >= 3)
    throw std::invalid_argument("dxreduce: cannot average all three axes (nothing to write)");

  // Determine remaining axis indices and labels.
  std::vector<int> rem_idx;
  for (int k = 0; k < 3; ++k)
    if (std::find(avg_idx.begin(), avg_idx.end(), k) == avg_idx.end())
      rem_idx.push_back(k);

  const char *latt_names[3] = {"a", "b", "c"};
  const char *cart_names[3] = {"x", "y", "z"};
  const bool cart = (atype == AxisType::Cartesian);

  // Build human-readable label strings.
  auto axis_label = [&](int k) { return cart ? cart_names[k] : latt_names[k]; };
  std::string avg_str, rem_str;
  for (int k : avg_idx) avg_str += axis_label(k);
  for (int k : rem_idx) rem_str += axis_label(k);

  trajan::log::info("dxreduce: averaging over {} → {}D output along {}",
                    avg_str, rem_idx.size(), rem_str);

  if (!cart) {
    // ── Lattice axes: fast exact path ────────────────────────────────────
    auto [data, remain] = average_lattice_axes(g, avg_idx);

    const occ::Vec3 *deltas[3] = {&g.da, &g.db, &g.dc};
    const size_t Ns[3] = {g.Na, g.Nb, g.Nc};

    if (remain.size() == 1) {
      const int r = remain[0];
      const double step = deltas[r]->norm();
      std::vector<double> pos(Ns[r]);
      for (size_t i = 0; i < Ns[r]; ++i) pos[i] = i * step;
      write_1d(opts.outfile, pos, data, latt_names[r], avg_str);
    } else {
      const int r0 = remain[0], r1 = remain[1];
      const double s0 = deltas[r0]->norm(), s1 = deltas[r1]->norm();
      std::vector<double> pos0(Ns[r0]), pos1(Ns[r1]);
      for (size_t i = 0; i < Ns[r0]; ++i) pos0[i] = i * s0;
      for (size_t i = 0; i < Ns[r1]; ++i) pos1[i] = i * s1;
      write_2d(opts.outfile, pos0, pos1, data,
               latt_names[r0], latt_names[r1], avg_str);
    }
  } else {
    // ── Cartesian axes ────────────────────────────────────────────────────
    if (is_orthogonal(g)) {
      // Orthogonal cell: x=a, y=b, z=c — use the fast lattice path.
      trajan::log::debug("dxreduce: orthogonal cell — Cartesian axes map directly to lattice axes");
      auto [data, remain] = average_lattice_axes(g, avg_idx);
      const occ::Vec3 *deltas[3] = {&g.da, &g.db, &g.dc};
      const size_t Ns[3] = {g.Na, g.Nb, g.Nc};
      if (remain.size() == 1) {
        const int r = remain[0];
        const double step = deltas[r]->norm();
        std::vector<double> pos(Ns[r]);
        for (size_t i = 0; i < Ns[r]; ++i) pos[i] = i * step;
        write_1d(opts.outfile, pos, data, cart_names[r], avg_str);
      } else {
        const int r0 = remain[0], r1 = remain[1];
        const double s0 = deltas[r0]->norm(), s1 = deltas[r1]->norm();
        std::vector<double> pos0(Ns[r0]), pos1(Ns[r1]);
        for (size_t i = 0; i < Ns[r0]; ++i) pos0[i] = i * s0;
        for (size_t i = 0; i < Ns[r1]; ++i) pos1[i] = i * s1;
        write_2d(opts.outfile, pos0, pos1, data,
                 cart_names[r0], cart_names[r1], avg_str);
      }
    } else {
      // Non-orthogonal cell: full Cartesian rebinning.
      trajan::log::debug("dxreduce: non-orthogonal cell — rebinning grid points by Cartesian coordinate");
      auto rb = rebin_cartesian(g, rem_idx);
      if (rem_idx.size() == 1) {
        write_1d(opts.outfile, rb.pos0, rb.values, cart_names[rem_idx[0]], avg_str);
      } else {
        write_2d(opts.outfile, rb.pos0, rb.pos1, rb.values,
                 cart_names[rem_idx[0]], cart_names[rem_idx[1]], avg_str);
      }
    }
  }

  trajan::log::info("dxreduce: wrote {}D output to {}", rem_idx.size(), opts.outfile);
}

CLI::App *add_dxreduce_subcommand(CLI::App &app) {
  CLI::App *cmd = app.add_subcommand(
      "dxreduce",
      "Load an OpenDX potential file and reduce its dimensionality by\n"
      "averaging over one or two axes. Writes a 1D or 2D potential profile.\n"
      "Supports both lattice axes (a,b,c) and Cartesian axes (x,y,z).");

  auto opts = std::make_shared<DXReduceOpts>();
  cmd->add_option("infile", opts->infile, "Input OpenDX file (.dx)")->required();
  cmd->add_option("-o,--output", opts->outfile, "Output file")->required();
  cmd->add_option(
         "--average", opts->average_axes,
         "Axes to collapse by averaging. Use lattice labels (a,b,c) or\n"
         "Cartesian labels (x,y,z) — do not mix the two sets.\n"
         "For orthogonal cells x=a, y=b, z=c (fast path).\n"
         "For non-orthogonal cells, Cartesian averaging uses rebinning.\n"
         "  --average c   → 2D map in the a-b plane\n"
         "  --average ab  → 1D profile along c\n"
         "  --average xy  → 1D profile along z\n"
         "  --average z   → 2D map in the x-y plane")
      ->required();

  cmd->callback([opts]() {
    trajan::log::set_subcommand_log_pattern("dxreduce");
    run_dxreduce_subcommand(*opts);
  });

  return cmd;
}

} // namespace trajan::main
