#pragma once
#include <omp.h>
#include <span>
#include <trajan/core/atom.h>
#include <trajan/core/linear_algebra.h>
#include <trajan/core/unit_cell.h>
#include <vector>

namespace trajan::core {

struct CellIndex {
  int a, b, c;
};

class Cell {
public:
  void add_atom(Atom &atom);
  const std::vector<Atom> &get_atoms() const { return m_atoms; }

private:
  std::vector<Atom> m_atoms;
  CellIndex m_index;
  friend class CellList;
};

struct CellListParameters {
  size_t num_ghosts;
  size_t a, b, c;
  size_t total_a, total_b, total_c;
  size_t a_upper, b_upper, c_upper;
  size_t a_end, b_end, c_end;
  CellListParameters(size_t a, size_t b, size_t c, size_t num_ghosts)
      : num_ghosts(num_ghosts), a(a), b(b), c(c), total_a(a + 2 * num_ghosts),
        total_b(b + 2 * num_ghosts), total_c(c + 2 * num_ghosts),
        a_end(a + num_ghosts), b_end(b + num_ghosts), c_end(c + num_ghosts),
        a_upper(a - num_ghosts), b_upper(b - num_ghosts),
        c_upper(c - num_ghosts) {}
};

class CellList {
public:
  CellList(const UnitCell &unit_cell, double cutoff, int num_threads);
  void update(std::span<Atom> atoms);
  template <typename Func> void for_each_pair(Func &&func) const;

private:
  static constexpr size_t NUMGHOSTS = 2;

  const UnitCell &m_unit_cell;
  const double m_cutoff, m_cutoffsq;
  const int m_num_threads;
  const CellListParameters m_params;
  std::vector<std::vector<std::vector<Cell>>> m_cells;
  const CellListParameters generate_cell_params(const UnitCell &unit_cell,
                                                double cutoff,
                                                const size_t num_ghosts) const;
  void initialise_cells();
  CellIndex get_cell_index(const Vec3 &pos) const;
  void clear_cells();
  template <typename Func>
  void process_cell(const Cell &center_cell, Func &&func) const;
  template <typename Func>
  void process_cell_pairs(const Cell &center_cell, int na, int nb, int nc,
                          Func &&func) const;
};

template <typename Func>
void CellList::process_cell_pairs(const Cell &center_cell, int na, int nb,
                                  int nc, Func &&func) const {
  auto center_atoms = center_cell.get_atoms();
  if (na == 0 && nb == 0 && nc == 0) {
    // for same cell, only check i < j pairs
    for (size_t i = 0; i < center_atoms.size(); ++i) {
      const auto &atom1 = center_atoms[i];
      for (size_t j = i + 1; j < center_atoms.size(); ++j) {
        const auto &atom2 = center_atoms[j];
        double rsq = atom1.square_distance(atom2);
        if (rsq <= m_cutoffsq) {
          func(atom1, atom2);
        }
      }
    }
    return;
  }

  // process forward direction to prevent double coutning
  if (na < 0 || (na == 0 && nb < 0) || (na == 0 && nb == 0 && nc < 0)) {
    return;
  }
  int neigh_a = center_cell.m_index.a + na;
  int neigh_b = center_cell.m_index.b + nb;
  int neigh_c = center_cell.m_index.c + nc;
  const Cell &neighbor_cell = m_cells[neigh_a][neigh_b][neigh_c];

  // process pairs between cells
  auto neigh_atoms = neighbor_cell.get_atoms();
  for (const Atom &atom1 : center_atoms) {
    for (const Atom &atom2 : neigh_atoms) {
      double rsq = atom1.square_distance(atom2);
      if (rsq <= m_cutoffsq) {
        func(atom1, atom2);
      }
    }
  }
}

template <typename Func>
void CellList::process_cell(const Cell &center_cell, Func &&func) const {
  const int adj = m_params.num_ghosts;
  for (int na = -adj; na <= adj; ++na) {
    for (int nb = -adj; nb <= adj; ++nb) {
      for (int nc = -adj; nc <= adj; ++nc) {
        process_cell_pairs(center_cell, na, nb, nc, func);
      }
    }
  }
}

template <typename Func> void CellList::for_each_pair(Func &&func) const {
  // #pragma omp parallel for collapse(3) schedule(dynamic) if (m_num_threads >
  // 1)
  for (int a = m_params.num_ghosts; a < m_params.a_end; ++a) {
    for (int b = m_params.num_ghosts; b < m_params.b_end; ++b) {
      for (int c = m_params.num_ghosts; c < m_params.c_end; ++c) {
        process_cell(m_cells[a][b][c], func);
      }
    }
  }
}
}; // namespace trajan::core

// #pragma once
// #include <memory>
// #include <omp.h>
// #include <span>
// #include <trajan/core/atom.h>
// #include <trajan/core/linear_algebra.h>
// #include <trajan/core/unit_cell.h>
// #include <unordered_map>
// #include <vector>
//
// namespace trajan::core {
//
// struct CellIndex {
//   size_t a, b, c;
// };
//
// class Cell {
// public:
//   void add_atom(Atom &atom);
//   const std::vector<Atom> &get_atoms() const { return m_atoms; }
//   void clear() { m_atoms.clear(); }
//
// private:
//   std::vector<Atom> m_atoms;
//   CellIndex m_index;
//   friend class CellList;
// };
//
// struct CellListParameters {
//   const size_t num_ghosts;
//   const size_t a, b, c;
//   const size_t total_cells;
//   const size_t total_a, total_b, total_c;
//   const size_t a_upper, b_upper, c_upper;
//   const size_t a_end, b_end, c_end;
//   CellListParameters(size_t a, size_t b, size_t c, size_t num_ghosts)
//       : num_ghosts(num_ghosts), a(a), b(b), c(c), total_a(a + 2 *
//       num_ghosts),
//         total_b(b + 2 * num_ghosts), total_c(c + 2 * num_ghosts),
//         total_cells(a + 2 * num_ghosts * b + 2 * num_ghosts * c +
//                     2 * num_ghosts),
//         a_end(a + num_ghosts), b_end(b + num_ghosts), c_end(c + num_ghosts),
//         a_upper(a - num_ghosts), b_upper(b - num_ghosts),
//         c_upper(c - num_ghosts) {}
// };
//
// struct CellIndexHash {
//   std::size_t operator()(const CellIndex &idx) const {
//     // Combine the three components into a single hash
//     // Using bit shifts to combine the values while minimizing collisions
//     return (static_cast<size_t>(idx.a) << 32) ^
//            (static_cast<size_t>(idx.b) << 16) ^ static_cast<size_t>(idx.c);
//   }
// };
//
// struct CellIndexEqual {
//   bool operator()(const CellIndex &lhs, const CellIndex &rhs) const {
//     return lhs.a == rhs.a && lhs.b == rhs.b && lhs.c == rhs.c;
//   }
// };
//
// class CellList {
// public:
//   static constexpr size_t NUMGHOSTS = 2;
//   static constexpr size_t NUMGHOSTSNEIGH =
//       ((2 * NUMGHOSTS + 1) * (2 * NUMGHOSTS + 1) * (2 * NUMGHOSTS + 1) - 1) /
//       2;
//   CellList(const UnitCell &unit_cell, double cutoff, int num_threads);
//   void update(std::span<Atom> atoms);
//   template <typename Func> void for_each_pair(Func &&func) const;
//
// private:
//   const UnitCell &m_unit_cell;
//   const double m_cutoff, m_cutoffsq;
//   const int m_num_threads;
//   const CellListParameters m_params;
//   std::vector<Cell> m_cells;
//   // std::unique_ptr<Cell[]> m_cells;
//
//   std::unordered_map<size_t, std::array<size_t, NUMGHOSTSNEIGH>>
//       m_neighbor_list;
//   std::unordered_map<size_t, size_t> m_neighbor_counts;
//   std::unordered_map<CellIndex, size_t, CellIndexHash, CellIndexEqual>
//       m_coordinate_map;
//   std::once_flag init_flag;
//
//   void initialise();
//   void build_coordinate_map();
//   void build_neighbor_list();
//   size_t get_index_from_coords(const CellIndex &coords) const;
//   constexpr size_t get_cell_index(int a, int b, int c) const {
//     return a * (m_params.total_b * m_params.total_c) + b * m_params.total_c +
//     c;
//   }
//   const CellListParameters generate_cell_params(const UnitCell &unit_cell,
//                                                 double cutoff,
//                                                 const size_t num_ghosts)
//                                                 const;
//   // void initialise_cells();
//   void clear_cells();
//   template <typename Func>
//   void process_cell(size_t cell_idx, Func &&func) const;
// };
// // template <typename Func>
// // void process_cell_pairs(const Cell &center_cell, int na, int nb, int nc,
// //                         Func &&func) const;
//
// // template <typename Func>
// // void CellList::process_cell_pairs(const Cell &center_cell, int na, int nb,
// //                                   int nc, Func &&func) const {
// //   auto center_atoms = center_cell.get_atoms();
// //   if (na == 0 && nb == 0 && nc == 0) {
// //     // for same cell, only check i < j pairs
// //     for (size_t i = 0; i < center_atoms.size(); ++i) {
// //       const auto &atom1 = center_atoms[i];
// //       for (size_t j = i + 1; j < center_atoms.size(); ++j) {
// //         const auto &atom2 = center_atoms[j];
// //         double rsq = atom1.square_distance(atom2);
// //         if (rsq <= m_cutoffsq) {
// //           func(atom1, atom2);
// //         }
// //       }
// //     }
// //     return;
// //   }
// //
// //   // process forward direction to prevent double coutning
// //   if (na < 0 || (na == 0 && nb < 0) || (na == 0 && nb == 0 && nc < 0)) {
// //     return;
// //   }
// //   // if (na < 0 || nb < 0 || nc < 0) {
// //   //   return;
// //   // }
// //   int neigh_a = center_cell.m_index.a + na;
// //   int neigh_b = center_cell.m_index.b + nb;
// //   int neigh_c = center_cell.m_index.c + nc;
// //   const Cell &neighbor_cell = m_cells[neigh_a][neigh_b][neigh_c];
// //
// //   // process pairs between cells
// //   auto neigh_atoms = neighbor_cell.get_atoms();
// //   for (const Atom &atom1 : center_atoms) {
// //     for (const Atom &atom2 : neigh_atoms) {
// //       double rsq = atom1.square_distance(atom2);
// //       if (rsq <= m_cutoffsq) {
// //         func(atom1, atom2);
// //       }
// //     }
// //   }
// // }
// //
// // template <typename Func>
// // void CellList::process_cell(const Cell &center_cell, Func &&func) const {
// //   const int adj = m_params.num_ghosts;
// //   for (int na = -adj; na <= adj; ++na) {
// //     for (int nb = -adj; nb <= adj; ++nb) {
// //       for (int nc = -adj; nc <= adj; ++nc) {
// //         process_cell_pairs(center_cell, na, nb, nc, func);
// //       }
// //     }
// //   }
// // }
// //
// // template <typename Func> void CellList::for_each_pair(Func &&func) const {
// //   // #pragma omp parallel for collapse(3) schedule(dynamic) if
// //   (m_num_threads
// //   >
// //   // 1)
// //   for (int a = m_params.num_ghosts; a < m_params.a_end; ++a) {
// //     for (int b = m_params.num_ghosts; b < m_params.b_end; ++b) {
// //       for (int c = m_params.num_ghosts; c < m_params.c_end; ++c) {
// //         process_cell(m_cells[a][b][c], func);
// //       }
// //     }
// //   }
// // }
// template <typename Func>
// void CellList::process_cell(size_t cell_idx, Func &&func) const {
//   const auto &center_cell = m_cells[cell_idx];
//   const auto &center_atoms = center_cell.get_atoms();
//
//   // Process atoms within the same cell
//   for (size_t i = 0; i < center_atoms.size(); ++i) {
//     const auto &atom1 = center_atoms[i];
//     for (size_t j = i + 1; j < center_atoms.size(); ++j) {
//       const auto &atom2 = center_atoms[j];
//       double rsq = atom1.square_distance(atom2);
//       if (rsq <= m_cutoffsq) {
//         func(atom1, atom2);
//       }
//     }
//   }
//
//   // Process neighboring cells
//   auto neighbor_it = m_neighbor_list.find(cell_idx);
//   if (neighbor_it != m_neighbor_list.end()) {
//     const auto &neighbor_indices = neighbor_it->second;
//     const size_t num_neighbors = m_neighbor_counts.at(cell_idx);
//
//     for (size_t i = 0; i < num_neighbors; ++i) {
//       const auto &neighbor_cell = m_cells[neighbor_indices[i]];
//       const auto &neighbor_atoms = neighbor_cell.get_atoms();
//
//       for (const auto &atom1 : center_atoms) {
//         for (const auto &atom2 : neighbor_atoms) {
//           double rsq = atom1.square_distance(atom2);
//           if (rsq <= m_cutoffsq) {
//             func(atom1, atom2);
//           }
//         }
//       }
//     }
//   }
// }
//
// template <typename Func> void CellList::for_each_pair(Func &&func) const {
//   // #pragma omp parallel for schedule(dynamic) if (m_num_threads > 1)
//   for (size_t idx = 0; idx < m_params.total_cells; ++idx) {
//     process_cell(idx, func);
//   }
// }
// }; // namespace trajan::core
