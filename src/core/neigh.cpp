#include <omp.h>
#include <stdexcept>
#include <trajan/core/neigh.h>
namespace trajan::core {

void Cell::add_atom(const Atom &atom) { m_atoms.push_back(atom); }

CellList::CellList(const UnitCell &unit_cell, double cutoff, int num_threads)
    : m_unit_cell(unit_cell), m_cutoff(cutoff), m_cutoffsq(cutoff * cutoff),
      m_num_threads(num_threads),
      m_params(generate_cell_params(unit_cell, cutoff, NUMGHOSTS)) {
  if (m_num_threads > 0) {
    omp_set_num_threads(m_num_threads);
  }
  {
#pragma omp critical
    initialise_cells();
  }
}

const CellListParameters
CellList::generate_cell_params(const UnitCell &unit_cell, double cutoff,
                               const size_t num_ghosts) const {
  return CellListParameters(
      static_cast<int>(
          std::floor(unit_cell.a_vector().norm() / (cutoff / num_ghosts))),
      static_cast<int>(
          std::floor(unit_cell.b_vector().norm() / (cutoff / num_ghosts))),
      static_cast<int>(
          std::floor(unit_cell.c_vector().norm() / (cutoff / num_ghosts))),
      num_ghosts);
}

size_t CellList::linear_index(size_t a, size_t b, size_t c) const {
  return (a * m_params.total_b * m_params.total_c) + (b * m_params.total_c) + c;
}

void CellList::initialise_cells() {
  m_cells.resize(m_params.total);
  for (size_t i = 0; i < m_params.total; i++) {
    const size_t a = i / (m_params.total_b * m_params.total_c);
    const size_t b =
        (i % (m_params.total_b * m_params.total_c)) / m_params.total_c;
    const size_t c = i % m_params.total_c;
    m_cells[i].m_index = CellIndex{a, b, c};
  }
  const int adj = static_cast<int>(m_params.num_ghosts);
  for (int a = m_params.num_ghosts; a < m_params.a_end; ++a) {
    for (int b = m_params.num_ghosts; b < m_params.b_end; ++b) {
      for (int c = m_params.num_ghosts; c < m_params.c_end; ++c) {
        std::vector<size_t> neighs;
        neighs.reserve(m_params.num_neighs);
        for (int na = -adj; na <= adj; ++na) {
          for (int nb = -adj; nb <= adj; ++nb) {
            for (int nc = -adj; nc <= adj; ++nc) {
              if (na < 0 || (na == 0 && nb < 0) ||
                  (na == 0 && nb == 0 && nc <= 0)) {
                continue;
              }
              neighs.push_back(
                  linear_index(static_cast<size_t>(static_cast<int>(a) + na),
                               static_cast<size_t>(static_cast<int>(b) + nb),
                               static_cast<size_t>(static_cast<int>(c) + nc)));
            }
          }
        }
        m_cell_neighs.insert(std::make_pair(linear_index(a, b, c), neighs));
      }
    }
  }
  m_cell_indices.reserve(m_params.total_real);
  for (const auto &pair : m_cell_neighs) {
    m_cell_indices.push_back(pair.first);
  }
}

void CellList::clear_cells() {
#pragma omp parallel for
  for (Cell &cell : m_cells) {
    cell.clear();
  }
}
void CellList::update(const std::span<Atom> atoms, const Mat3N atoms_pos) {
  if (atoms.size() != atoms_pos.cols()) {
    std::runtime_error("Number of atoms and their positions not the same.");
  }
  clear_cells();
  Mat3N frac_pos = m_unit_cell.to_fractional(atoms_pos);
  frac_pos = frac_pos.array() - frac_pos.array().floor();
  IVec inds_a = (frac_pos.row(0) * m_params.a).cast<int>();
  IVec inds_b = (frac_pos.row(1) * m_params.b).cast<int>();
  IVec inds_c = (frac_pos.row(2) * m_params.c).cast<int>();
  Mat3N pos = m_unit_cell.to_cartesian(frac_pos);
  for (int atom_i = 0; atom_i < atoms.size(); atom_i++) {
    Atom atom = atoms[atom_i];
    Vec3 atom_pos = pos.col(atom_i);
    atom.update_position(atom_pos);
    int ind_a = inds_a[atom_i];
    int ind_b = inds_b[atom_i];
    int ind_c = inds_c[atom_i];
    cell_at(ind_a + m_params.num_ghosts, ind_b + m_params.num_ghosts,
            ind_c + m_params.num_ghosts)
        .add_atom(atom);
    if (ind_a >= m_params.num_ghosts && ind_a < m_params.a_upper &&
        ind_b >= m_params.num_ghosts && ind_b < m_params.b_upper &&
        ind_c >= m_params.num_ghosts && ind_c < m_params.c_upper) {
      continue;
    }
    int a = (ind_a < m_params.num_ghosts)
                ? 1
                : (ind_a >= m_params.a_upper ? -1 : 0);
    int b = (ind_b < m_params.num_ghosts)
                ? 1
                : (ind_b >= m_params.b_upper ? -1 : 0);
    int c = (ind_c < m_params.num_ghosts)
                ? 1
                : (ind_c >= m_params.c_upper ? -1 : 0);

    for (int ia = 0; ia <= std::abs(a); ia++) {
      for (int ib = 0; ib <= std::abs(b); ib++) {
        for (int ic = 0; ic <= std::abs(c); ic++) {
          if (ia == 0 && ib == 0 && ic == 0) {
            continue;
          }
          int a_shift = ia * a, b_shift = ib * b, c_shift = ic * c;
          Vec3 shift = m_unit_cell.direct() * Vec3(a_shift, b_shift, c_shift);
          cell_at(ind_a + a_shift * m_params.a + m_params.num_ghosts,
                  ind_b + b_shift * m_params.b + m_params.num_ghosts,
                  ind_c + c_shift * m_params.c + m_params.num_ghosts)
              .m_atoms.emplace_back(atom.create_ghost(shift));
        }
      }
    }
  }
}
} // namespace trajan::core
