#include <optional>
#include <stdexcept>
#include <trajan/core/graph.h>
#include <trajan/core/log.h>
#include <trajan/core/molecule.h>
#include <trajan/core/neigh.h>
#include <trajan/core/trajectory.h>
#include <trajan/core/unit_cell.h>
#include <trajan/io/file_handler.h>
#include <trajan/io/selection.h>
#include <vector>

namespace trajan::core {

Trajectory::~Trajectory() {
  for (auto &handler : m_handlers) {
    if (handler) {
      handler->finalise();
    }
  }
}

void Trajectory::load_files(const std::vector<fs::path> &files) {
  m_handlers = io::read_input_files(files);

  if (m_handlers.empty()) {
    throw std::runtime_error("Could not load files.");
  }

  m_handlers[0]->initialise();
}

bool Trajectory::next_frame() {
  if (m_handlers.empty()) {
    return false;
  }

  if (m_current_handler_index >= m_handlers.size()) {
    throw std::runtime_error("Incorrect index.");
  }
  if (m_handlers[m_current_handler_index]->read_frame(m_frame)) {
    m_frame_loaded = true;
    m_current_frame_index++;
    return true;
  }

  m_handlers[m_current_handler_index]->finalise();
  m_current_handler_index++;

  while (m_current_handler_index < m_handlers.size()) {
    if (m_handlers[m_current_handler_index]->initialise()) {
      if (m_handlers[m_current_handler_index]->read_frame(m_frame)) {
        m_frame_loaded = true;
        m_current_frame_index++;
        return true;
      }
    }
    m_handlers[m_current_handler_index]->finalise();
    m_current_handler_index++;
  }
  m_frame_loaded = false;

  return false;
}

void Trajectory::reset() {
  for (auto &handler : m_handlers) {
    if (handler) {
      handler->finalise();
    }
  }

  m_current_handler_index = 0;
  m_current_frame_index = 0;
  m_frame_loaded = false;
}

void Trajectory::update_neigh_uc(const UnitCell &unit_cell) {
  // TODO: compare with current uc to see if it needs updating
  m_topo_neigh_list.update_uc(unit_cell);
}

void Trajectory::update_neigh_rcut(double rcut) {
  // TODO: compare with current rcut to see if it needs updating
  m_topo_neigh_list.update_rcut(rcut);
}

void Trajectory::update_neigh_threads(size_t threads) {
  // TODO: compare with current threads to see if it needs updating
  m_topo_neigh_list.update_threads(threads);
}

void Trajectory::update_neigh() {
  m_topo_neigh_list.update_uc(this->unit_cell());
}

void Trajectory::update_neigh(const UnitCell &unit_cell, double rcut,
                              size_t threads) {
  m_topo_neigh_list = NeighbourList(unit_cell, rcut, threads);
}

std::vector<EntityType>
Trajectory::get_entities(const io::SelectionCriteria &selection) {
  if (!m_frame_loaded) {
    throw std::runtime_error("No frame loaded. Call next_frame() first.");
  }
  std::vector<Atom> atoms = this->atoms();
  std::vector<Molecule> molecules;
  std::vector<EntityType> entities;
  entities.reserve(atoms.size());

  trajan::log::debug("Processing selection {}", selection.index());

  if (std::holds_alternative<io::MoleculeIndexSelection>(selection) ||
      std::holds_alternative<io::MoleculeTypeSelection>(selection)) {
    // build molecules
    molecules = this->unit_cell_molecules();
  }

  entities = std::visit(
      [&](const auto &sel) {
        using ActualType = std::decay_t<decltype(sel)>;
        return io::process_selection<ActualType>(sel, atoms, molecules,
                                                 entities);
      },
      selection);

  size_t num_entities = entities.size();
  if (num_entities == 0) {
    // FIXME: Fix the selection names for clarity.
    throw std::runtime_error("No entities found in selection.");
  }
  trajan::log::debug("Identified {} entities in selection", num_entities);

  return entities;
}

std::vector<EntityType>
Trajectory::get_entities(const std::vector<io::SelectionCriteria> &selections) {
  std::vector<EntityType> entities;
  for (const io::SelectionCriteria &sel : selections) {
    std::vector<EntityType> e2 = this->get_entities(sel);
    entities.reserve(entities.size() + e2.size());
    entities.insert(entities.end(), e2.begin(), e2.end());
  }
  return entities;
}

void Trajectory::update_bond_graph() {
  std::vector<Atom> atoms = this->atoms();
  m_bond_graph = BondGraph(atoms);
}

void Trajectory::update_bond_graph(std::vector<Atom> &atoms) {
  m_bond_graph = BondGraph(atoms);
}

const BondGraph Trajectory::unit_cell_connectivity() {
  if (m_unit_cell_connectivity_needs_update) {
    update_unit_cell_connectivity();
  }
  return m_bond_graph;
}

void Trajectory::update_unit_cell_connectivity() {
  this->update_neigh();
  this->update_bond_graph();

  m_topo_neigh_list.update(this->atoms());

  // TODO: allow user input tolerance
  double tol = 0.4;
  const std::vector<Atom> &atoms = this->atoms();
  NeighbourCallback func = [&](const Entity &ent1, const Entity &ent2,
                               double rsq) {
    const Atom &atom1 = atoms[ent1.idx];
    const Atom &atom2 = atoms[ent2.idx];
    // std::visit(
    //     [rsq, &bg, tol](const auto &e1, const auto &e2) {
    //       using T1 = std::decay_t<decltype(e1)>;
    //       using T2 = std::decay_t<decltype(e2)>;
    //       if constexpr (std::is_same_v<T1, Atom> && std::is_same_v<T2, Atom>)
    //       {
    //         std::optional<Bond> bond = e1.is_bonded_with_rsq(e2, rsq, tol);
    //         if (bond) {
    //           bg.add_edge(e1.index, e2.index, *bond);
    //         }
    //       }
    //     },
    //     ntt1, ntt2);
    std::optional<Bond> bond = atom1.is_bonded_with_rsq(atom2, rsq, tol);
    if (bond) {
      m_bond_graph.add_edge(atom1.index, atom2.index, *bond);
    }
  };

  this->m_topo_neigh_list.iterate_neighbours(func);

  m_unit_cell_connectivity_needs_update = false;
}

void Trajectory::update_unit_cell_molecules() {
  using CC = graph::ConnectedComponent<Atom, Bond>;
  m_unit_cell_molecules.clear();
  BondGraph bg = unit_cell_connectivity();
  std::vector<CC> ccs = bg.find_connected_components();
  for (CC &cc : ccs) {
    m_unit_cell_molecules.emplace_back(cc);
  }
  m_unit_cell_molecules_needs_update = false;
}

const std::vector<Molecule> &Trajectory::unit_cell_molecules() {
  if (m_unit_cell_molecules_needs_update) {
    update_unit_cell_molecules();
  }
  return m_unit_cell_molecules;
}

} // namespace trajan::core
