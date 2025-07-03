#pragma once
#include <filesystem>
#include <trajan/core/frame.h>
#include <trajan/core/graph.h>
#include <trajan/core/neigh.h>
#include <trajan/core/topology.h>
#include <trajan/core/unit_cell.h>
#include <trajan/io/file_handler.h>
#include <trajan/io/selection.h>
#include <vector>

namespace trajan::core {

namespace fs = std::filesystem;

class Trajectory {

public:
  Trajectory(const Trajectory &) = delete;
  Trajectory &operator=(const Trajectory &) = delete;

  Trajectory(Trajectory &&) = default;
  Trajectory &operator=(Trajectory &&) = default;

  Trajectory() = default;
  ~Trajectory();

  void load_files(const std::vector<fs::path> &files);

  bool next_frame();

  void reset();

  inline bool has_frames() const {
    return !m_handlers.empty() && m_current_handler_index < m_handlers.size();
  }

  inline size_t current_frame_index() const { return m_current_frame_index; }

  inline Frame &frame() { return m_frame; }

  inline const std::vector<Atom> &atoms() const { return m_frame.atoms(); }

  inline size_t num_atoms() const { return m_frame.num_atoms(); }

  inline const UnitCell &unit_cell() const { return m_frame.unit_cell(); }

  const std::vector<Molecule> extract_molecules();

  std::vector<EntityType> get_entities(const io::SelectionCriteria &selection);

  std::vector<EntityType>
  get_entities(const std::vector<io::SelectionCriteria> &selections);

  const Topology &get_topology();

  void update_topology();

private:
  bool m_guess_connectivity{true};

  size_t m_current_handler_index{0};
  size_t m_current_frame_index{0};
  bool m_frame_loaded{false};
  std::vector<io::FileHandlerPtr> m_handlers;
  Frame m_frame;

  bool advance_to_next_available_frame();

  // void update_neigh();
  // void update_neigh(const UnitCell &unit_cell, double rcut, size_t threads);
  // void update_neigh_uc(const UnitCell &unit_cell);
  // void update_neigh_rcut(double rcut);
  // void update_neigh_threads(size_t threads);
  // void update_unit_cell_topology();
  // void update_unit_cell_molecules();

  mutable Topology m_topology;
  mutable std::vector<Molecule> m_molecules{};
  mutable bool m_topology_needs_update{true};
};

}; // namespace trajan::core
