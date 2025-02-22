// src file
#include <omp.h>
#include <trajan/core/atom.h>
#include <trajan/core/neigh.h>
#include <trajan/core/unit_cell.h>
#include <vector>

namespace trajan::core {

void Cell::add_atom(Atom &atom) {
  std::lock_guard<std::mutex> lock(*m_mutex);
  m_atoms.emplace_back(atom);
}

CellList::CellList(const UnitCell &unit_cell, double cutoff, int num_threads)
    : m_unit_cell(unit_cell), m_cutoff(cutoff), m_cell_size(cutoff / 2.0),
      m_ghost_cells(3), m_num_threads(num_threads) {
  if (m_num_threads > 0) {
    omp_set_num_threads(m_num_threads);
  }
  initialize_cells();
}

void CellList::initialize_cell_indices() {
  for (int x = 0; x < m_cells.size(); ++x) {
    for (int y = 0; y < m_cells[0].size(); ++y) {
      for (int z = 0; z < m_cells[0][0].size(); ++z) {
        m_cells[x][y][z].m_index = CellIndex{x, y, z};
      }
    }
  }
}

void CellList::initialize_cells() {
  auto [nx, ny, nz] = calculate_cell_counts();
  const int total_x = nx + 2 * m_ghost_cells;
  const int total_y = ny + 2 * m_ghost_cells;
  const int total_z = nz + 2 * m_ghost_cells;

  m_cells.resize(total_x);
  for (auto &plane : m_cells) {
    plane.resize(total_y);
    for (auto &line : plane) {
      line.resize(total_z);
    }
  }

  initialize_cell_indices();
}

// void CellList::update(std::span<Atom> atoms) {
//   clear_cells();
//   // std::for_each(std::execution::par_unseq, atoms.begin(), atoms.end(),
//   //               [this](Atom &atom) { assign_atom_to_cell(atom); });
// #pragma omp parallel for schedule(dynamic)
//   for (size_t i = 0; i < atoms.size(); ++i) {
//     assign_atom_to_cell(atoms[i]);
//   }
// }

void CellList::update(std::span<Atom> atoms) {
  clear_cells();

  // Pre-calculate cell assignments to avoid lock contention
  struct AtomCell {
    Atom *atom;
    CellIndex idx;
  };
  std::vector<AtomCell> assignments;
  assignments.reserve(atoms.size());

#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < atoms.size(); ++i) {
    auto idx = get_cell_indices(atoms[i].position());
#pragma omp critical
    assignments.push_back({&atoms[i], idx});
  }

  // Bulk insert atoms into cells
  for (const auto &assignment : assignments) {
    m_cells[assignment.idx.x][assignment.idx.y][assignment.idx.z].add_atom(
        *assignment.atom);
    handle_periodic_images(*assignment.atom);
  }
}

void CellList::clear_cells() {
  const size_t nx = m_cells.size();
  const size_t ny = m_cells[0].size();
  const size_t nz = m_cells[0][0].size();

// Use collapse(3) to parallelize all three nested loops
#pragma omp parallel for collapse(3) schedule(static)
  for (size_t x = 0; x < nx; ++x) {
    for (size_t y = 0; y < ny; ++y) {
      for (size_t z = 0; z < nz; ++z) {
        std::lock_guard<std::mutex> lock(*m_cells[x][y][z].m_mutex);
        m_cells[x][y][z].m_atoms.clear();
      }
    }
  }
}

CellIndex CellList::get_cell_indices(const Eigen::Vector3d &position) const {
  Eigen::Vector3d frac_pos = m_unit_cell.to_fractional(position);
  return CellIndex{static_cast<int>(frac_pos.x() / m_cell_size) + m_ghost_cells,
                   static_cast<int>(frac_pos.y() / m_cell_size) + m_ghost_cells,
                   static_cast<int>(frac_pos.z() / m_cell_size) +
                       m_ghost_cells};
}

void CellList::assign_atom_to_cell(Atom &atom) {
  auto indices = get_cell_indices(atom.position());
  m_cells[indices.x][indices.y][indices.z].add_atom(atom);
  handle_periodic_images(atom);
}

// void CellList::handle_periodic_images(Atom &atom) {
//   for (int x = -1; x <= 1; ++x) {
//     for (int y = -1; y <= 1; ++y) {
//       for (int z = -1; z <= 1; ++z) {
//         if (x == 0 && y == 0 && z == 0)
//           continue;
//         create_periodic_image(atom, x, y, z);
//       }
//     }
//   }
// }

void CellList::handle_periodic_images(Atom &atom) {
  Vec3 pos = atom.position();
  Vec3 frac_pos = m_unit_cell.to_fractional(pos);

  // Only create periodic images near boundaries
  for (int x = -1; x <= 1; ++x) {
    if (x != 0 && std::abs(frac_pos.x() - 0.5) < m_cell_size)
      continue;
    for (int y = -1; y <= 1; ++y) {
      if (y != 0 && std::abs(frac_pos.y() - 0.5) < m_cell_size)
        continue;
      for (int z = -1; z <= 1; ++z) {
        if (z != 0 && std::abs(frac_pos.z() - 0.5) < m_cell_size)
          continue;
        if (x == 0 && y == 0 && z == 0)
          continue;
        create_periodic_image(atom, x, y, z);
      }
    }
  }
}

void CellList::create_periodic_image(Atom &atom, int x, int y, int z) {
  Eigen::Vector3d shift = m_unit_cell.direct() * Eigen::Vector3d(x, y, z);
  Eigen::Vector3d image_pos = atom.position() + shift;

  auto indices = get_cell_indices(image_pos);
  if (is_ghost_cell(indices)) {
    m_cells[indices.x][indices.y][indices.z].add_atom(atom);
  }
}

bool CellList::is_ghost_cell(const CellIndex &indices) const {
  return indices.x < m_ghost_cells ||
         indices.x >= m_cells.size() - m_ghost_cells ||
         indices.y < m_ghost_cells ||
         indices.y >= m_cells[0].size() - m_ghost_cells ||
         indices.z < m_ghost_cells ||
         indices.z >= m_cells[0][0].size() - m_ghost_cells;
}

double CellList::calculate_minimum_image_distance(const Atom &atom1,
                                                  const Atom &atom2) const {
  Eigen::Vector3d dr = atom2.position() - atom1.position();
  Eigen::Vector3d frac_dr = m_unit_cell.to_fractional(dr);
  frac_dr = frac_dr.array() - (frac_dr.array().round());
  dr = m_unit_cell.to_cartesian(frac_dr);
  return dr.norm();
}

std::tuple<int, int, int> CellList::calculate_cell_counts() const {
  return {
      static_cast<int>(std::ceil(m_unit_cell.a_vector().norm() / m_cell_size)),
      static_cast<int>(std::ceil(m_unit_cell.b_vector().norm() / m_cell_size)),
      static_cast<int>(std::ceil(m_unit_cell.c_vector().norm() / m_cell_size))};
}

} // namespace trajan::core
