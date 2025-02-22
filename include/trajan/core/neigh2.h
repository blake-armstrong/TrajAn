// header file
#pragma once

#include <compare>
#include <mutex>
#include <span>
#include <trajan/core/atom.h>
#include <trajan/core/linear_algebra.h>
#include <trajan/core/unit_cell.h>
#include <vector>

namespace trajan::core {

struct CellIndex {
  int x, y, z;

  bool operator==(const CellIndex &) const = default;
  auto operator<=>(const CellIndex &) const = default;
};

class Cell {
public:
  Cell() : m_mutex(std::make_unique<std::mutex>()) {}

  // Explicitly define move constructor and assignment operator
  Cell(Cell &&other) noexcept
      : m_atoms(std::move(other.m_atoms)), m_mutex(std::move(other.m_mutex)),
        m_index(other.m_index) {}

  Cell &operator=(Cell &&other) noexcept {
    if (this != &other) {
      m_atoms = std::move(other.m_atoms);
      m_mutex = std::move(other.m_mutex);
      m_index = other.m_index;
    }
    return *this;
  }

  // Delete copy operations
  Cell(const Cell &) = delete;
  Cell &operator=(const Cell &) = delete;

  void add_atom(Atom &atom);
  const std::vector<std::reference_wrapper<Atom>> &get_atoms() const {
    return m_atoms;
  }

private:
  std::vector<std::reference_wrapper<Atom>> m_atoms;
  std::unique_ptr<std::mutex> m_mutex; // Wrap mutex in unique_ptr
  CellIndex m_index;
  friend class CellList; // Allows CellList to access private members
};

class CellList {
public:
  CellList(const UnitCell &unit_cell, double cutoff, int num_threads);

  // Updates the cell lists with new atom positions
  void update(std::span<Atom> atoms);

  // Executes a function for all pairs of atoms within cutoff distance
  template <typename Func> void for_each_pair(Func &&func) const;

private:
  const UnitCell &m_unit_cell;
  const double m_cutoff;
  const double m_cell_size;
  const int m_ghost_cells;
  std::vector<std::vector<std::vector<Cell>>> m_cells;
  const int m_num_threads;

  void initialize_cells();
  void initialize_cell_indices();
  void clear_cells();
  CellIndex get_cell_indices(const Vec3 &position) const;
  void assign_atom_to_cell(Atom &atom);
  void handle_periodic_images(Atom &atom);
  void create_periodic_image(Atom &atom, int x, int y, int z);
  bool is_ghost_cell(const CellIndex &indices) const;
  std::tuple<int, int, int> calculate_cell_counts() const;

  template <typename Func>
  void process_neighbor_cell(const Cell &center_cell,
                             const CellIndex &center_idx, int dx, int dy,
                             int dz, Func &&func) const;
  template <typename Func>
  void process_cell_pairs(const CellIndex &cell_idx, Func &&func) const;

  double calculate_minimum_image_distance(const Atom &atom1,
                                          const Atom &atom2) const;
};

// template <typename Func>
// void CellList::process_neighbor_cell(const Cell &center_cell,
//                                      const CellIndex &center_idx, int dx,
//                                      int dy, int dz, Func &&func) const {
//   // Calculate the neighbor cell indices with periodic wrapping
//   // First convert sizes to signed integers to avoid signed/unsigned mismatch
//   const int nx = static_cast<int>(m_cells.size());
//   const int ny = static_cast<int>(m_cells[0].size());
//   const int nz = static_cast<int>(m_cells[0][0].size());
//
//   // Calculate neighbor indices with periodic wrapping
//   // Use intermediate calculations to ensure proper handling of negative
//   values CellIndex neigh_idx{((center_idx.x + dx + nx) % nx + nx) % nx,
//                       ((center_idx.y + dy + ny) % ny + ny) % ny,
//                       ((center_idx.z + dz + nz) % nz + nz) % nz};
//
//   const Cell &neighbor_cell = m_cells[neigh_idx.x][neigh_idx.y][neigh_idx.z];
//
//   // Skip empty cells
//   if (center_cell.get_atoms().empty() || neighbor_cell.get_atoms().empty()) {
//     return;
//   }
//
//   // Process pairs between center cell and neighbor cell
//   for (const auto &atom1_ref : center_cell.get_atoms()) {
//     const Atom &atom1 = atom1_ref.get();
//
//     for (const auto &atom2_ref : neighbor_cell.get_atoms()) {
//       const Atom &atom2 = atom2_ref.get();
//
//       // Skip self-interaction when processing the same cell
//       if (&atom1 == &atom2) {
//         continue;
//       }
//
//       // Calculate the minimum image distance between atoms
//       double distance = calculate_minimum_image_distance(atom1, atom2);
//
//       // If atoms are within cutoff, apply the user-provided function
//       if (distance <= m_cutoff) {
//         func(atom1, atom2);
//       }
//     }
//   }
// }

template <typename Func>
void CellList::process_neighbor_cell(const Cell &center_cell,
                                     const CellIndex &center_idx, int dx,
                                     int dy, int dz, Func &&func) const {
  const int nx = static_cast<int>(m_cells.size());
  const int ny = static_cast<int>(m_cells[0].size());
  const int nz = static_cast<int>(m_cells[0][0].size());

  CellIndex neigh_idx{((center_idx.x + dx + nx) % nx + nx) % nx,
                      ((center_idx.y + dy + ny) % ny + ny) % ny,
                      ((center_idx.z + dz + nz) % nz + nz) % nz};

  // Skip if we would double count
  if (dx == 0 && dy == 0 && dz == 0) {
    // For same cell, only check i < j pairs
    for (size_t i = 0; i < center_cell.get_atoms().size(); ++i) {
      const auto &atom1 = center_cell.get_atoms()[i].get();
      for (size_t j = i + 1; j < center_cell.get_atoms().size(); ++j) {
        const auto &atom2 = center_cell.get_atoms()[j].get();
        double distance = calculate_minimum_image_distance(atom1, atom2);
        if (distance <= m_cutoff) {
          func(atom1, atom2);
        }
      }
    }
    return;
  }

  // For different cells, only process if we're in the forward direction
  // This prevents double counting pairs between cells
  if (dx < 0 || (dx == 0 && dy < 0) || (dx == 0 && dy == 0 && dz < 0)) {
    return;
  }

  const Cell &neighbor_cell = m_cells[neigh_idx.x][neigh_idx.y][neigh_idx.z];

  // Skip empty cells
  if (center_cell.get_atoms().empty() || neighbor_cell.get_atoms().empty()) {
    return;
  }

  // Process pairs between cells
  for (const auto &atom1_ref : center_cell.get_atoms()) {
    const Atom &atom1 = atom1_ref.get();
    for (const auto &atom2_ref : neighbor_cell.get_atoms()) {
      const Atom &atom2 = atom2_ref.get();
      double distance = calculate_minimum_image_distance(atom1, atom2);
      if (distance <= m_cutoff) {
        func(atom1, atom2);
      }
    }
  }
}

template <typename Func>
void CellList::process_cell_pairs(const CellIndex &cell_idx,
                                  Func &&func) const {
  const auto &center_cell = m_cells[cell_idx.x][cell_idx.y][cell_idx.z];
  const int range = static_cast<int>(std::ceil(m_cutoff / m_cell_size));

  for (int dx = -range; dx <= range; ++dx) {
    for (int dy = -range; dy <= range; ++dy) {
      for (int dz = -range; dz <= range; ++dz) {
        process_neighbor_cell(center_cell, cell_idx, dx, dy, dz,
                              std::forward<Func>(func));
      }
    }
  }
}

template <typename Func> void CellList::for_each_pair(Func &&func) const {
  const int real_start = m_ghost_cells;
  const int x_end = m_cells.size() - m_ghost_cells;
  const int y_end = m_cells[0].size() - m_ghost_cells;
  const int z_end = m_cells[0][0].size() - m_ghost_cells;

#pragma omp parallel for collapse(3) schedule(dynamic) if (m_num_threads > 1)

  for (int x = real_start; x < x_end; ++x) {
    for (int y = real_start; y < y_end; ++y) {
      for (int z = real_start; z < z_end; ++z) {
        process_cell_pairs(CellIndex{x, y, z}, func);
      }
    }
  }
}

} // namespace trajan::core
