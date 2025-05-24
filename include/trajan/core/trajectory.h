#pragma once
#include <filesystem>
#include <trajan/core/frame.h>
#include <trajan/core/graph.h>
#include <trajan/core/neigh.h>
#include <trajan/core/unit_cell.h>
#include <trajan/io/file_handler.h>
#include <trajan/io/selection.h>
#include <vector>

namespace trajan::core {

namespace fs = std::filesystem;

class Trajectory {

public:
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

  void update_bond_graph();
  void update_bond_graph(std::vector<Atom> &atoms);
  const std::vector<Molecule> &unit_cell_molecules();

  std::vector<EntityType> get_entities(const io::SelectionCriteria &selection);
  std::vector<EntityType>
  get_entities(const std::vector<io::SelectionCriteria> &selections);

private:
  bool m_guess_connectivity{true};

  size_t m_current_handler_index{0};
  size_t m_current_frame_index{0};
  bool m_frame_loaded{false};
  std::vector<io::FileHandlerPtr> m_handlers;
  Frame m_frame;

  bool advance_to_next_available_frame();

  void update_neigh();
  void update_neigh(const UnitCell &unit_cell, double rcut, size_t threads);
  void update_neigh_uc(const UnitCell &unit_cell);
  void update_neigh_rcut(double rcut);
  void update_neigh_threads(size_t threads);

  const BondGraph unit_cell_connectivity();
  void update_unit_cell_connectivity();
  void update_unit_cell_molecules();

  NeighbourList m_topo_neigh_list;
  BondGraph m_bond_graph;
  std::vector<BondGraph::NodeID> m_bond_graph_ids;
  std::vector<Molecule> m_unit_cell_molecules{};
  bool m_unit_cell_connectivity_needs_update{true};
  bool m_unit_cell_molecules_needs_update{true};
};

}; // namespace trajan::core
