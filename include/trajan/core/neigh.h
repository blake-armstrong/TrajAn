#pragma once
#include <ankerl/unordered_dense.h>
#include <fmt/format.h>
#include <omp.h>
#include <stdexcept>
#include <trajan/core/atom.h>
#include <trajan/core/linear_algebra.h>
#include <trajan/core/log.h>
#include <trajan/core/molecule.h>
#include <trajan/core/unit_cell.h>
#include <variant>
#include <vector>

namespace trajan::core {

struct Entity {
  size_t idx;
  double x, y, z;
  enum class Type { Atom, Molecule } type;

  Entity(size_t idx, const Vec3 &pos, Type t)
      : idx(idx), x(pos.x()), y(pos.y()), z(pos.z()), type(t) {}
  Entity(const Atom &atom)
      : idx(atom.index), x(atom.x), y(atom.y), z(atom.z), type(Type::Atom) {}
  Entity(const Molecule &molecule)
      : idx(molecule.index), x(molecule.x), y(molecule.y), z(molecule.z),
        type(Type::Molecule) {}

  inline Vec3 position() const { return {x, y, z}; }
  inline double square_distance(const Entity &other) const {
    double dx = other.x - x, dy = other.y - y, dz = other.z - z;
    return dx * dx + dy * dy + dz * dz;
  }

  bool operator==(const Atom &atom) const { return type == Type::Atom; }
  bool operator==(const Molecule &molecule) const {
    return type == Type::Molecule;
  }
};

using EntityType = std::variant<Atom, Molecule>;

struct VariantHash {
  std::size_t operator()(const std::variant<Atom, Molecule> &var) const {
    if (std::holds_alternative<Atom>(var)) {
      const Atom &atom = std::get<Atom>(var);
      return std::hash<int>{}(atom.index);
    } else {
      const Molecule &molecule = std::get<Molecule>(var);
      return std::hash<int>{}(molecule.index);
    }
  }
};

struct VariantEqual {
  bool operator()(const std::variant<Atom, Molecule> &lhs,
                  const std::variant<Atom, Molecule> &rhs) const {
    if (lhs.index() != rhs.index()) {
      return false;
    }

    if (std::holds_alternative<Atom>(lhs)) {
      return std::get<Atom>(lhs) == std::get<Atom>(rhs);
    } else {
      return std::get<Molecule>(lhs) == std::get<Molecule>(rhs);
    }
  }
};

struct NeighbourListPacket {
  Mat3N wrapped_cart_pos;
  Mat3N wrapped_frac_pos;
  std::vector<Entity::Type> obj_types;
  std::vector<size_t> presence_tracker;
  bool check_presence{false};
  ankerl::unordered_dense::map<size_t, size_t> index_to_canonical;

  NeighbourListPacket() = default;
  NeighbourListPacket(const std::vector<Entity::Type> &obj_types,
                      const Mat3N &cart_pos, const Mat3N &frac_pos) {
    this->wrapped_cart_pos = cart_pos;
    this->wrapped_frac_pos = frac_pos;
    this->obj_types = obj_types;
    this->check();
  }

  size_t size() const { return obj_types.size(); }

  void check() const {
    if (this->wrapped_cart_pos.size() != this->wrapped_frac_pos.size()) {
      throw std::runtime_error(
          "Cartesian and fractional positions must be the same size");
    };
    if (this->wrapped_cart_pos.cols() != this->obj_types.size()) {
      throw std::runtime_error(
          "Entities vector and position matrix not same size.");
    };
  }

  inline bool are_same_entity(size_t idx1, size_t idx2) const {
    auto it1 = index_to_canonical.find(idx1);
    auto it2 = index_to_canonical.find(idx2);

    if (it1 == index_to_canonical.end() || it2 == index_to_canonical.end()) {
      return false;
    }

    return it1->second == it2->second;
  }
};

using NeighbourCallback =
    std::function<void(const Entity &, const Entity &, double)>;

// base class with common interface
class NeighbourListBase {
public:
  // NeighbourListBase(UnitCell &unit_cell, double rcut, size_t num_threads);
  virtual ~NeighbourListBase() = default;

  NeighbourListBase(const UnitCell &unit_cell, double cutoff,
                    size_t num_threads);
  virtual void update(const NeighbourListPacket &np) = 0;
  virtual void iterate_neighbours(const NeighbourCallback &callback) const = 0;

  virtual void update_uc(const UnitCell &unit_cell) = 0;
  virtual void update_rcut(double rcut) = 0;
  virtual void update_threads(size_t threads) = 0;

  inline const UnitCell &unit_cell() { return m_unit_cell; }

protected:
  const UnitCell &m_unit_cell;
  double m_cutoff, m_cutoffsq;
  size_t m_num_threads;
};

struct CellIndex {
  size_t a, b, c;
};

class CellList; // forward declaration

class Cell {
public:
  inline void add_entity(const Entity &entity) { m_entities.push_back(entity); }
  inline void add_entity(const size_t idx, Vec3 &pos, Entity::Type type) {
    m_entities.emplace_back(idx, pos, type);
  }
  const std::vector<Entity> &get_entities() const { return m_entities; }
  void clear() { m_entities.clear(); }

private:
  std::vector<Entity> m_entities;
  CellIndex m_index;
  friend class CellList;
};

struct CellListParameters {
  size_t num_ghosts;
  size_t num_neighs;
  size_t a, b, c;
  size_t total_a, total_b, total_c;
  size_t total;
  size_t total_real;
  size_t a_upper, b_upper, c_upper;
  size_t a_end, b_end, c_end;
  CellListParameters() = default;
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

// cell list algorithm for neighbours
class CellList : public NeighbourListBase {
public:
  CellList(const UnitCell &unit_cell, double cutoff, size_t num_threads);
  void update(const NeighbourListPacket &np) override;
  void iterate_neighbours(const NeighbourCallback &callback) const override {
    this->cell_loop(callback);
  };

  void update_uc(const UnitCell &unit_cell) override {
    CellList(unit_cell, m_cutoff, m_num_threads);
  }
  void update_rcut(double rcut) override {
    CellList(m_unit_cell, rcut, m_num_threads);
  }
  void update_threads(size_t num_threads) override {
    CellList(m_unit_cell, m_cutoff, num_threads);
  }

private:
  // TODO: test non-periodic
  static constexpr size_t CELLDIVISOR = 2;
  static constexpr size_t GHOSTCELLS = CELLDIVISOR;

  NeighbourListPacket m_current_nlp;
  bool m_dummy{false};
  CellListParameters m_params;
  std::vector<Cell> m_cells;
  ankerl::unordered_dense::map<size_t, std::vector<size_t>> m_cell_neighs;
  std::vector<size_t> m_cell_indices;

  inline size_t linear_index(size_t a, size_t b, size_t c) const {
    return (a * m_params.total_b * m_params.total_c) + (b * m_params.total_c) +
           c;
  }
  inline Cell &cell_at(size_t a, size_t b, size_t c) {
    return m_cells[linear_index(a, b, c)];
  }
  inline const Cell &cell_at(size_t a, size_t b, size_t c) const {
    return m_cells[linear_index(a, b, c)];
  }

  CellListParameters
  generate_cell_params(size_t ghost_cells = GHOSTCELLS) const;
  size_t determine_opt_num_threads();
  void initialise_cells(size_t ghost_cells = GHOSTCELLS);
  void clear_cells();
  void cell_loop(const NeighbourCallback &callback) const;
};

// verlet/double loop algorithm for neighbours
class VerletList : public NeighbourListBase {
public:
  VerletList(const UnitCell &unit_cell, double cutoff, size_t num_threads);
  void update(const NeighbourListPacket &np) override;
  void iterate_neighbours(const NeighbourCallback &callback) const override {
    this->verlet_loop(callback);
  };

  void update_uc(const UnitCell &unit_cell) override {
    VerletList(unit_cell, m_cutoff, m_num_threads);
  }
  void update_rcut(double rcut) override {
    VerletList(m_unit_cell, rcut, m_num_threads);
  }
  void update_threads(size_t num_threads) override {
    VerletList(m_unit_cell, m_cutoff, num_threads);
  }

private:
  NeighbourListPacket m_current_nlp;
  size_t determine_opt_num_threads();
  void verlet_loop(const NeighbourCallback &callback) const;
};

// runtime-selectable NeighbourList class
class NeighbourList {
public:
  enum class Type { Cell, Verlet };

  NeighbourList() = default;
  NeighbourList(const UnitCell &unit_cell, double cutoff,
                size_t num_threads = 1, Type type = Type::Cell) {

    switch (type) {
    case Type::Cell:
      m_impl = std::make_unique<CellList>(unit_cell, cutoff, num_threads);
      break;
    case Type::Verlet:
      m_impl = std::make_unique<VerletList>(unit_cell, cutoff, num_threads);
      break;
    }
    m_init = true;
  }

  void update(const std::vector<Atom> &atoms);
  void update(const std::vector<EntityType> &og_objects,
              Molecule::Origin o = Molecule::CentreOfMass);
  void update(const std::vector<std::vector<EntityType>> &og_objects_vec,
              Molecule::Origin o = Molecule::CentreOfMass);

  void iterate_neighbours(const NeighbourCallback &callback) {
    m_impl->iterate_neighbours(callback);
  }

  void update_uc(const UnitCell &unit_cell) { m_impl->update_uc(unit_cell); }
  void update_rcut(double rcut) { m_impl->update_rcut(rcut); }
  void update_threads(size_t num_threads) {
    m_impl->update_threads(num_threads);
  }

private:
  std::vector<EntityType> m_original_objects;
  NeighbourListPacket m_current_nlp;
  std::unique_ptr<NeighbourListBase> m_impl;
  bool m_init{false};

  void base_update(const std::vector<EntityType> &og_objects,
                   Molecule::Origin o);
};

}; // namespace trajan::core
