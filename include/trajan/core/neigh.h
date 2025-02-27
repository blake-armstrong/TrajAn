#pragma once
#include <omp.h>
#include <span>
#include <trajan/core/atom.h>
#include <trajan/core/linear_algebra.h>
#include <trajan/core/unit_cell.h>
#include <unordered_map>
#include <vector>

namespace trajan::core {

// cell list algorithm for neighbours
struct CellIndex {
  size_t a, b, c;
};

class Cell {
public:
  void add_atom(const Atom &atom);
  const std::vector<Atom> &get_atoms() const { return m_atoms; }
  void clear() { m_atoms.clear(); }

private:
  std::vector<Atom> m_atoms;
  CellIndex m_index;
  friend class CellList;
};

struct CellListParameters {
  const size_t num_ghosts;
  const size_t num_neighs;
  const size_t a, b, c;
  const size_t total_a, total_b, total_c;
  const size_t total;
  const size_t total_real;
  const size_t a_upper, b_upper, c_upper;
  const size_t a_end, b_end, c_end;
  CellListParameters(size_t a, size_t b, size_t c, size_t num_ghosts)
      : num_ghosts(num_ghosts),
        num_neighs(((2 * num_ghosts + 1) * (2 * num_ghosts + 1) *
                        (2 * num_ghosts + 1) -
                    1) /
                   2),
        a(a), b(b), c(c), total_a(a + 2 * num_ghosts),
        total_b(b + 2 * num_ghosts), total_c(c + 2 * num_ghosts),
        total((a + 2 * num_ghosts) * (b + 2 * num_ghosts) *
              (c + 2 * num_ghosts)),
        total_real(a * b * c), a_end(a + num_ghosts), b_end(b + num_ghosts),
        c_end(c + num_ghosts), a_upper(a - num_ghosts), b_upper(b - num_ghosts),
        c_upper(c - num_ghosts) {}
};

class CellList {
public:
  CellList(const UnitCell &unit_cell, double cutoff, int num_threads);
  void update(const std::span<Atom> atoms, Mat3N atoms_pos);
  template <typename Func> void cell_loop(Func &&func) const;

private:
  static constexpr size_t NUMGHOSTS = 2;

  const UnitCell &m_unit_cell;
  const double m_cutoff, m_cutoffsq;
  const int m_num_threads;
  const CellListParameters m_params;
  std::vector<Cell> m_cells;
  std::unordered_map<size_t, std::vector<size_t>> m_cell_neighs;
  std::vector<size_t> m_cell_indices;

  const CellListParameters generate_cell_params(const UnitCell &unit_cell,
                                                double cutoff,
                                                const size_t num_ghosts) const;
  size_t linear_index(size_t a, size_t b, size_t c) const;
  Cell &cell_at(size_t a, size_t b, size_t c) {
    return m_cells[linear_index(a, b, c)];
  }
  const Cell &cell_at(size_t a, size_t b, size_t c) const {
    return m_cells[linear_index(a, b, c)];
  }
  void initialise_cells();
  void clear_cells();
};

template <typename Func> void CellList::cell_loop(Func &&func) const {
#pragma omp parallel for num_threads(m_num_threads)
  for (size_t cell_i = 0; cell_i < m_params.total_real; cell_i++) {
    size_t center_cell_idx = m_cell_indices[cell_i];
    auto center_atoms = m_cells[center_cell_idx].get_atoms();
    for (size_t i = 0; i < center_atoms.size(); ++i) {
      const Atom &atom1 = center_atoms[i];
      for (size_t j = i + 1; j < center_atoms.size(); ++j) {
        const Atom &atom2 = center_atoms[j];
        double rsq = atom1.square_distance(atom2);
        if (rsq <= m_cutoffsq) {
          func(atom1, atom2, rsq);
        }
      }
    }
    std::vector<size_t> neighs = m_cell_neighs.at(center_cell_idx);
    for (const size_t &neigh_idx : neighs) {
      auto neigh_atoms = m_cells[neigh_idx].get_atoms();
      for (const Atom &atom1 : center_atoms) {
        for (const Atom &atom2 : neigh_atoms) {
          double rsq = atom1.square_distance(atom2);
          if (rsq <= m_cutoffsq) {
            func(atom1, atom2, rsq);
          }
        }
      }
    }
  }
}

// verlet/double loop algorithm for neighbours
class VerletList {
public:
  VerletList(const UnitCell &unit_cell, double cutoff, int num_threads);
  void update(const std::span<Atom> atoms, Mat3N atoms_pos);
  template <typename Func> void verlet_loop(Func &&func) const;

private:
  const UnitCell &m_unit_cell;
  const double m_cutoff, m_cutoffsq;
  const int m_num_threads;
  std::span<Atom> m_atoms;
  Mat3N m_atom_positions;
};

template <typename Func> void VerletList::verlet_loop(Func &&func) const {
#pragma omp parallel for collapse(2) num_threads(m_num_threads)
  for (size_t i = 0; i < m_atoms.size(); ++i) {
    for (size_t j = i + 1; j < m_atoms.size(); ++j) {
      Vec3 dr = m_atoms[j].position() - m_atoms[i].position();
      Vec3 frac_dr = m_unit_cell.to_fractional(dr);
      frac_dr = frac_dr.array() - frac_dr.array().round();
      dr = m_unit_cell.to_cartesian(frac_dr);
      double rsq = dr.squaredNorm();
      if (rsq <= m_cutoffsq) {
        func(m_atoms[i], m_atoms[j], rsq);
      }
    }
  }
};

// generic neighbour class:
// allows easy switching between algorithms
using NeighbourCallback =
    std::function<void(const Atom &, const Atom &, double)>;

class INeighbourList {
public:
  virtual ~INeighbourList() = default;
  virtual void update(const std::span<Atom> atoms, Mat3N atoms_pos) = 0;
  virtual void iterate_neighbours(const NeighbourCallback &callback) const = 0;
};

class CellListNeighbours : public INeighbourList {
public:
  CellListNeighbours(const UnitCell &unit_cell, double cutoff, int num_threads)
      : m_cell_list(unit_cell, cutoff, num_threads) {}
  void update(const std::span<Atom> atoms, Mat3N atoms_pos) override {
    m_cell_list.update(atoms, atoms_pos);
  }
  void iterate_neighbours(const NeighbourCallback &callback) const override {
    m_cell_list.cell_loop(callback);
  }

private:
  CellList m_cell_list;
};

class VerletListNeighbours : public INeighbourList {
public:
  VerletListNeighbours(const UnitCell &unit_cell, double cutoff,
                       int num_threads)
      : m_verlet_list(unit_cell, cutoff, num_threads) {}
  void update(const std::span<Atom> atoms, Mat3N atoms_pos) override {
    m_verlet_list.update(atoms, atoms_pos);
  }
  void iterate_neighbours(const NeighbourCallback &callback) const override {
    m_verlet_list.verlet_loop(callback);
  }

private:
  VerletList m_verlet_list;
};

enum class NeighbourListType { CELL_LIST, VERLET_LIST };

class NeighbourListHandler {
public:
  NeighbourListHandler(const UnitCell &unit_cell, double cutoff,
                       int num_threads,
                       NeighbourListType type = NeighbourListType::CELL_LIST)
      : m_unit_cell(unit_cell), m_cutoff(cutoff), m_num_threads(num_threads) {
    set_list_type(type);
  }
  void set_list_type(NeighbourListType type) {
    switch (type) {
    case NeighbourListType::CELL_LIST:
      m_neighbour_list = std::make_unique<CellListNeighbours>(
          m_unit_cell, m_cutoff, m_num_threads);
      break;
    case NeighbourListType::VERLET_LIST:
      m_neighbour_list = std::make_unique<VerletListNeighbours>(
          m_unit_cell, m_cutoff, m_num_threads);
      break;
    }
  }
  void set_list_type(const std::string &type_str) {
    if (type_str == "cell" || type_str == "cell_list") {
      set_list_type(NeighbourListType::CELL_LIST);
    } else if (type_str == "verlet" || type_str == "verlet_list") {
      set_list_type(NeighbourListType::VERLET_LIST);
    } else {
      throw std::runtime_error("Unknown neighbour list type: " + type_str);
    }
  }
  void update(const std::span<Atom> atoms, Mat3N atoms_pos) {
    m_neighbour_list->update(atoms, atoms_pos);
  }
  void iterate_neighbours(const NeighbourCallback &callback) const {
    m_neighbour_list->iterate_neighbours(callback);
  }
  template <typename Func> void for_each_neighbour(Func &&func) const {
    NeighbourCallback callback = std::forward<Func>(func);
    iterate_neighbours(callback);
  }

private:
  const UnitCell &m_unit_cell;
  double m_cutoff;
  int m_num_threads;
  std::unique_ptr<INeighbourList> m_neighbour_list;
};

}; // namespace trajan::core
