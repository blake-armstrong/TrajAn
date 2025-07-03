#include <optional>
#include <stdexcept>
#include <trajan/core/element.h>
#include <trajan/core/graph.h>
#include <trajan/core/log.h>
#include <trajan/core/molecule.h>
#include <trajan/core/neigh.h>
#include <trajan/core/topology.h>
#include <trajan/core/trajectory.h>
#include <trajan/core/unit_cell.h>
#include <trajan/io/file_handler.h>
#include <trajan/io/selection.h>
#include <vector>

namespace trajan::core {

using runtime_values::max_cov_cutoff;

Trajectory::~Trajectory() {
  for (auto &handler : m_handlers) {
    if (handler) {
      trajan::log::critical("hello");
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

// void Trajectory::update_neigh_uc(const UnitCell &unit_cell) {
//   // TODO: compare with current uc to see if it needs updating
//   m_topo_neigh_list.update_uc(unit_cell);
// }
//
// void Trajectory::update_neigh_rcut(double rcut) {
//   // TODO: compare with current rcut to see if it needs updating
//   m_topo_neigh_list.update_rcut(rcut);
// }
//
// void Trajectory::update_neigh_threads(size_t threads) {
//   // TODO: compare with current threads to see if it needs updating
//   m_topo_neigh_list.update_threads(threads);
// }
//
// void Trajectory::update_neigh() {
//   m_topo_neigh_list.update_uc(this->unit_cell());
// }
//
// void Trajectory::update_neigh(const UnitCell &unit_cell, double rcut,
//                               size_t threads) {
//   m_topo_neigh_list = NeighbourList(unit_cell, rcut, threads);
// }

std::vector<EntityType>
Trajectory::get_entities(const io::SelectionCriteria &selection) {
  std::vector<Atom> atoms = this->atoms();
  std::vector<Molecule> molecules;
  std::vector<EntityType> entities;
  entities.reserve(atoms.size());

  trajan::log::debug("Processing selection {}", selection.index());

  if (std::holds_alternative<io::MoleculeIndexSelection>(selection) ||
      std::holds_alternative<io::MoleculeTypeSelection>(selection)) {
    // build molecules
    molecules = this->extract_molecules();
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

const Topology &Trajectory::get_topology() {
  if (m_topology_needs_update) {
    this->update_topology();
  }
  return m_topology;
}

void Trajectory::update_topology() {
  const std::vector<Atom> &atoms = this->atoms();
  trajan::log::debug("update_topology: Starting with {} atoms", atoms.size());

  if (atoms.empty()) {
    trajan::log::error("No atoms available for topology update");
    m_topology_needs_update = false;
    return;
  }

  trajan::log::debug("Creating BondGraph...");
  BondGraph bg(atoms);
  trajan::log::debug("BondGraph created successfully");
  size_t threads = 1; // FIXME:
  trajan::log::debug("Creating NeighbourList with cutoff {:.3f}",
                     max_cov_cutoff);
  NeighbourList nl(this->unit_cell(), max_cov_cutoff, threads);
  trajan::log::debug("NeighbourList created, updating with atoms...");

  nl.update(atoms);
  trajan::log::debug("NeighbourList updated successfully");

  double tol = 0.4; // FIXME:
  size_t bond_count = 0;

  NeighbourCallback func = [&](const Entity &ent1, const Entity &ent2,
                               double rsq) {
    const Atom &atom1 = atoms[ent1.idx];
    const Atom &atom2 = atoms[ent2.idx];
    std::optional<Bond> bond = atom1.is_bonded_with_rsq(atom2, rsq, tol);
    if (bond) {
      bg.add_edge(atom1.index, atom2.index, bond.value());
      bond_count++;
    }
  };

  trajan::log::debug("Starting neighbour iteration...");
  nl.iterate_neighbours(func);
  trajan::log::debug("Neighbour iteration complete, found {} bonds",
                     bond_count);
  m_topology = Topology(bg);
}

const std::vector<Molecule> Trajectory::extract_molecules() {
  if (m_topology_needs_update) {
    this->update_topology();
  }
  return m_topology.extract_molecules();
}

} // namespace trajan::core
