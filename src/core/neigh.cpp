#include <omp.h>
#include <trajan/core/neigh.h>
namespace trajan::core {

void Cell::add_atom(Atom &atom) { m_atoms.push_back(atom); }

CellList::CellList(const UnitCell &unit_cell, double cutoff, int num_threads)
    : m_unit_cell(unit_cell), m_cutoff(cutoff), m_cutoffsq(cutoff * cutoff),
      m_num_threads(num_threads),
      m_params(generate_cell_params(unit_cell, cutoff, NUMGHOSTS)) {
  if (m_num_threads > 0) {
    omp_set_num_threads(m_num_threads);
  }
  initialise_cells();
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

void CellList::initialise_cells() {
  m_cells.resize(m_params.total_a);
  for (auto &plane : m_cells) {
    plane.resize(m_params.total_b);
    for (auto &line : plane) {
      line.resize(m_params.total_c);
    }
  }
  for (int a = 0; a < m_params.total_a; ++a) {
    for (int b = 0; b < m_params.total_b; ++b) {
      for (int c = 0; c < m_params.total_c; ++c) {
        m_cells[a][b][c].m_index = CellIndex{a, b, c};
      }
    }
  }
}

void CellList::clear_cells() {
  // #pragma omp parallel for collapse(3) schedule(static)
  for (size_t a = 0; a < m_params.total_a; ++a) {
    for (size_t b = 0; b < m_params.total_b; ++b) {
      for (size_t c = 0; c < m_params.total_c; ++c) {
        m_cells[a][b][c].m_atoms.clear();
      }
    }
  }
}
void CellList::update(std::span<Atom> atoms) {
  clear_cells();
  // #pragma omp parallel for schedule(dynamic)
  for (Atom &atom : atoms) {
    Vec3 frac_pos = m_unit_cell.to_fractional(atom.position());
    // frac_pos = frac_pos.array() - frac_pos.array().floor();
    int ind_a = static_cast<int>(frac_pos.x() * m_params.a);
    int ind_b = static_cast<int>(frac_pos.y() * m_params.b);
    int ind_c = static_cast<int>(frac_pos.z() * m_params.c);
    m_cells[ind_a + m_params.num_ghosts][ind_b + m_params.num_ghosts]
           [ind_c + m_params.num_ghosts]
               .add_atom(atom);

    // do the ghosts as well
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
          m_cells[ind_a + a_shift * m_params.a + m_params.num_ghosts]
                 [ind_b + b_shift * m_params.b + m_params.num_ghosts]
                 [ind_c + c_shift * m_params.c + m_params.num_ghosts]
                     .m_atoms.emplace_back(atom.create_ghost(shift));
        }
      }
    }
  }
}
} // namespace trajan::core

// #include <iostream>
// #include <memory>
// #include <omp.h>
// #include <trajan/core/neigh.h>
// namespace trajan::core {
//
// void Cell::add_atom(Atom &atom) { m_atoms.push_back(atom); }
//
// CellList::CellList(const UnitCell &unit_cell, double cutoff, int num_threads)
//     : m_unit_cell(unit_cell), m_cutoff(cutoff), m_cutoffsq(cutoff * cutoff),
//       m_num_threads(num_threads),
//       m_params(generate_cell_params(unit_cell, cutoff, NUMGHOSTS)) {
//   // Ensure thread-safe initialization
//   std::call_once(init_flag, [this]() { initialise(); });
// }
//
// void CellList::initialise() {
//   // Set thread count first
//   if (m_num_threads > 0) {
//     omp_set_num_threads(m_num_threads);
//   }
//
//   // Pre-allocate memory for all containers
//   m_cells = std::vector<Cell>(m_params.total_cells);
//   m_coordinate_map.reserve(m_params.total_cells);
//   m_neighbor_list.reserve(m_params.total_cells);
//   m_neighbor_counts.reserve(m_params.total_cells);
//
//   // Build coordinate map without parallelization
//   build_coordinate_map();
//
//   // Build neighbor list without parallelization
//   build_neighbor_list();
// }
//
// void CellList::build_coordinate_map() {
//   // for (size_t idx = 0; idx < m_params.total_cells; ++idx) {
//   // }
//   for (size_t a = 0; a < m_params.total_a; ++a) {
//     for (size_t b = 0; b < m_params.total_b; ++b) {
//       for (size_t c = 0; c < m_params.total_c; ++c) {
//         CellIndex idx{a, b, c};
//         size_t linear_idx = get_cell_index(a, b, c);
//         m_coordinate_map[idx] = linear_idx;
//         m_cells[linear_idx].m_index = idx;
//       }
//     }
//   }
// }
//
// void CellList::build_neighbor_list() {
//   const int adj = m_params.num_ghosts;
//
//   for (int a = m_params.num_ghosts; a < m_params.a_end; ++a) {
//     for (int b = m_params.num_ghosts; b < m_params.b_end; ++b) {
//       for (int c = m_params.num_ghosts; c < m_params.c_end; ++c) {
//         size_t center_idx = get_cell_index(a, b, c);
//         size_t neighbor_count = 0;
//         auto &neighbors = m_neighbor_list[center_idx];
//
//         for (int na = -adj; na <= adj; ++na) {
//           for (int nb = -adj; nb <= adj; ++nb) {
//             for (int nc = -adj; nc <= adj; ++nc) {
//               // Skip self and backward directions to prevent double counting
//               if (na < 0 || (na == 0 && nb < 0) ||
//                   (na == 0 && nb == 0 && nc <= 0)) {
//                 continue;
//               }
//
//               int neigh_a = a + na;
//               int neigh_b = b + nb;
//               int neigh_c = c + nc;
//
//               size_t neighbor_idx = get_cell_index(neigh_a, neigh_b,
//               neigh_c); neighbors[neighbor_count++] = neighbor_idx;
//             }
//           }
//         }
//
//         m_neighbor_counts[center_idx] = neighbor_count;
//       }
//     }
//   }
// }
//
// size_t CellList::get_index_from_coords(const CellIndex &coords) const {
//   auto it = m_coordinate_map.find(coords);
//   if (it != m_coordinate_map.end()) {
//     return it->second;
//   }
//   // Handle error case - could throw exception or return sentinel value
//   throw std::out_of_range("Invalid cell coordinates");
// }
//
// const CellListParameters
// CellList::generate_cell_params(const UnitCell &unit_cell, double cutoff,
//                                const size_t num_ghosts) const {
//   return CellListParameters(
//       static_cast<int>(
//           std::floor(unit_cell.a_vector().norm() / (cutoff / num_ghosts))),
//       static_cast<int>(
//           std::floor(unit_cell.b_vector().norm() / (cutoff / num_ghosts))),
//       static_cast<int>(
//           std::floor(unit_cell.c_vector().norm() / (cutoff / num_ghosts))),
//       num_ghosts);
// }
//
// // void CellList::initialise_cells() {
// //   // m_params = calculate_cell_counts(NUMGHOSTS);
// //   m_cells.resize(m_params.total_a);
// //   for (auto &plane : m_cells) {
// //     plane.resize(m_params.total_b);
// //     for (auto &line : plane) {
// //       line.resize(m_params.total_c);
// //     }
// //   }
// //   for (int a = 0; a < m_params.total_a; ++a) {
// //     for (int b = 0; b < m_params.total_b; ++b) {
// //       for (int c = 0; c < m_params.total_c; ++c) {
// //         m_cells[a][b][c].m_index = CellIndex{a, b, c};
// //       }
// //     }
// //   }
// // }
//
// // void CellList::clear_cells() {
// //   // #pragma omp parallel for collapse(3) schedule(static)
// //   for (size_t a = 0; a < m_params.total_a; ++a) {
// //     for (size_t b = 0; b < m_params.total_b; ++b) {
// //       for (size_t c = 0; c < m_params.total_c; ++c) {
// //         m_cells[a][b][c].m_atoms.clear();
// //       }
// //     }
// //   }
// // }
//
// void CellList::clear_cells() {
//   // #pragma omp parallel for schedule(static) if (m_num_threads > 1)
//   for (size_t idx = 0; idx < m_params.total_cells; ++idx) {
//     m_cells[idx].clear();
//   }
// }
//
// // void CellList::update(std::span<Atom> atoms) {
// //   clear_cells();
// //   // #pragma omp parallel for schedule(dynamic)
// //   for (Atom &atom : atoms) {
// //     Vec3 frac_pos = m_unit_cell.to_fractional(atom.position());
// //     // frac_pos = frac_pos.array() - frac_pos.array().floor();
// //     int ind_a = static_cast<int>(frac_pos.x() * m_params.a);
// //     int ind_b = static_cast<int>(frac_pos.y() * m_params.b);
// //     int ind_c = static_cast<int>(frac_pos.z() * m_params.c);
// //     m_cells[ind_a + m_params.num_ghosts][ind_b + m_params.num_ghosts]
// //            [ind_c + m_params.num_ghosts]
// //                .add_atom(atom);
// //
// //     // do the ghosts as well
// //     if (ind_a >= m_params.num_ghosts && ind_a < m_params.a_upper &&
// //         ind_b >= m_params.num_ghosts && ind_b < m_params.b_upper &&
// //         ind_c >= m_params.num_ghosts && ind_c < m_params.c_upper) {
// //       continue;
// //     }
// //     int a = (ind_a < m_params.num_ghosts)
// //                 ? 1
// //                 : (ind_a >= m_params.a_upper ? -1 : 0);
// //     int b = (ind_b < m_params.num_ghosts)
// //                 ? 1
// //                 : (ind_b >= m_params.b_upper ? -1 : 0);
// //     int c = (ind_c < m_params.num_ghosts)
// //                 ? 1
// //                 : (ind_c >= m_params.c_upper ? -1 : 0);
// //
// //     for (int ia = 0; ia <= std::abs(a); ia++) {
// //       for (int ib = 0; ib <= std::abs(b); ib++) {
// //         for (int ic = 0; ic <= std::abs(c); ic++) {
// //           if (ia == 0 && ib == 0 && ic == 0) {
// //             continue;
// //           }
// //           int a_shift = ia * a, b_shift = ib * b, c_shift = ic * c;
// //           Vec3 shift = m_unit_cell.direct() * Vec3(a_shift, b_shift,
// //           c_shift); m_cells[ind_a + a_shift * m_params.a +
// //           m_params.num_ghosts]
// //                  [ind_b + b_shift * m_params.b + m_params.num_ghosts]
// //                  [ind_c + c_shift * m_params.c + m_params.num_ghosts]
// //                      .m_atoms.emplace_back(atom.create_ghost(shift));
// //         }
// //       }
// //     }
// //   }
// // }
// void CellList::update(std::span<Atom> atoms) {
//   clear_cells();
//
//   for (Atom &atom : atoms) {
//     Vec3 frac_pos = m_unit_cell.to_fractional(atom.position());
//     size_t ind_a = static_cast<int>(frac_pos.x() * m_params.a);
//     size_t ind_b = static_cast<int>(frac_pos.y() * m_params.b);
//     size_t ind_c = static_cast<int>(frac_pos.z() * m_params.c);
//
//     // Use the coordinate map to get the cell index
//     CellIndex coords{ind_a + m_params.num_ghosts, ind_b +
//     m_params.num_ghosts,
//                      ind_c + m_params.num_ghosts};
//     std::cout << "aa" << std::endl;
//     size_t cell_idx;
//     {
// #pragma omp critical
//       cell_idx = get_index_from_coords(coords);
//       // std::cout << "bb: " << cell_idx << std::endl;
//       // std::cout << "m_cells: " << m_params.total_cells << std::endl;
//       m_cells[cell_idx].add_atom(atom);
//     }
//     std::cout << "cc" << std::endl;
//
//     // Handle ghost cells...
//     if (ind_a >= m_params.num_ghosts && ind_a < m_params.a_upper &&
//         ind_b >= m_params.num_ghosts && ind_b < m_params.b_upper &&
//         ind_c >= m_params.num_ghosts && ind_c < m_params.c_upper) {
//       continue;
//     }
//
//     int a = (ind_a < m_params.num_ghosts)
//                 ? 1
//                 : (ind_a >= m_params.a_upper ? -1 : 0);
//     int b = (ind_b < m_params.num_ghosts)
//                 ? 1
//                 : (ind_b >= m_params.b_upper ? -1 : 0);
//     int c = (ind_c < m_params.num_ghosts)
//                 ? 1
//                 : (ind_c >= m_params.c_upper ? -1 : 0);
//
//     for (int ia = 0; ia <= std::abs(a); ia++) {
//       for (int ib = 0; ib <= std::abs(b); ib++) {
//         for (int ic = 0; ic <= std::abs(c); ic++) {
//           if (ia == 0 && ib == 0 && ic == 0) {
//             continue;
//           }
//           int a_shift = ia * a, b_shift = ib * b, c_shift = ic * c;
//           Vec3 shift = m_unit_cell.direct() * Vec3(a_shift, b_shift,
//           c_shift);
//
//           CellIndex ghost_coords{
//               ind_a + a_shift * static_cast<int>(m_params.a) +
//                   static_cast<int>(m_params.num_ghosts),
//               ind_b + b_shift * static_cast<int>(m_params.b) +
//                   static_cast<int>(m_params.num_ghosts),
//               ind_c + c_shift * static_cast<int>(m_params.c) +
//                   static_cast<int>(m_params.num_ghosts)};
//
//           size_t ghost_idx = get_index_from_coords(ghost_coords);
//           m_cells[ghost_idx].m_atoms.emplace_back(atom.create_ghost(shift));
//         }
//       }
//     }
//   }
// }
// } // namespace trajan::core
