#pragma once
#include <trajan/core/frame.h>
#include <trajan/core/graph.h>
#include <trajan/core/neigh.h>
#include <trajan/core/unit_cell.h>
#include <trajan/io/selection.h>
#include <vector>

namespace trajan::core {

class Trajectory {

public:
  Trajectory() = default;
  // Trajectory(Frame &frame);

  inline Frame &frame() { return m_frame; }
  inline const std::vector<Atom> &atoms() const { return m_frame.atoms(); }
  inline size_t num_atoms() const { return m_frame.num_atoms(); }

  inline const UnitCell &unit_cell() const { return m_frame.unit_cell(); }
  // void update_frame();

  void update_neigh();
  void update_neigh(const UnitCell &unit_cell, double rcut, size_t threads);
  void update_neigh_uc(const UnitCell &unit_cell);
  void update_neigh_rcut(double rcut);
  void update_neigh_threads(size_t threads);

  void update_bond_graph();
  void update_bond_graph(std::vector<Atom> &atoms);
  const std::vector<Molecule> &unit_cell_molecules();

  Entities get_entities();
  Entities get_entities(const io::SelectionCriteria &selection);
  // Entities get_entities(const std::vector<io::SelectionCriteria>
  // &selections);
  NeighbourListPacket
  get_neighpack_from_entities(const Entities &entities,
                              Molecule::Origin o = Molecule::CentreOfMass);

private:
  // NeighbourList m_neigh_list;
  bool m_guess_connectivity{true};

  const BondGraph unit_cell_connectivity();
  void update_unit_cell_connectivity();
  void update_unit_cell_molecules();

  NeighbourList m_topo_neigh_list;
  BondGraph m_bond_graph;
  Frame m_frame;
  std::vector<BondGraph::NodeID> m_bond_graph_ids;
  std::vector<Molecule> m_unit_cell_molecules{};
  bool m_unit_cell_connectivity_needs_update{true};
  bool m_unit_cell_molecules_needs_update{true};
};

}; // namespace trajan::core
