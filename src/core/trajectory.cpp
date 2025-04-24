#include "trajan/core/graph.h"
#include "trajan/core/molecule.h"
#include <stdexcept>
#include <trajan/core/log.h>
#include <trajan/core/neigh.h>
#include <trajan/core/trajectory.h>
#include <trajan/io/selection.h>
#include <vector>

namespace trajan::core {

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

Entities Trajectory::get_entities() {
  std::vector<Atom> atoms = this->atoms();
  Entities entities;
  entities.reserve(atoms.size());
  for (Atom &atom : atoms) {
    entities.push_back(atom);
  }
  return entities;
}

Entities Trajectory::get_entities(const io::SelectionCriteria &selection) {
  using SelType = std::decay_t<decltype(selection)>;
  std::vector<Atom> atoms = this->atoms();
  std::vector<Molecule> molecules;
  Entities entities;
  entities.reserve(atoms.size());

  trajan::log::debug("Processing selection {}", typeid(SelType).name());

  if constexpr (std::is_same_v<SelType, io::MoleculeSelection>) {
    // build molecules
    molecules = this->unit_cell_molecules();
  }
  entities =
      io::process_selection<SelType>(selection, atoms, molecules, entities);

  return entities;
}

// Entities
// Trajectory::get_entities(const std::vector<io::SelectionCriteria>
// &selections) {
//   Entities entities;
//   for (const io::SelectionCriteria &sel : selections) {
//     Entities e2 = this->get_entities(sel);
//     entities.reserve(entities.size() + e2.size());
//     entities.insert(entities.end(), e2.begin(), e2.end());
//   }
//   return entities;
// }

NeighbourListPacket
Trajectory::get_neighpack_from_entities(const Entities &entities,
                                        Molecule::Origin o) {
  NeighbourListPacket nlp;
  bool all_atoms =
      std::all_of(entities.begin(), entities.end(), [](const auto &entity) {
        return std::holds_alternative<Atom>(entity);
      });
  if (all_atoms) {
    Frame frame = this->frame();
    return NeighbourListPacket(entities, frame.wrapped_cart_pos(),
                               frame.frac_pos());
  }

  size_t e_size = entities.size();
  Mat3N cart_pos(3, e_size);
  for (size_t i = 0; i < e_size; i++) {
    const EntityType &ent = entities[i];
    std::visit(
        [&cart_pos, &o, i](const auto &entity) {
          using T = std::decay_t<decltype(entity)>;
          if constexpr (std::is_same_v<T, Atom>) {
            cart_pos(0, i) = entity.x;
            cart_pos(1, i) = entity.y;
            cart_pos(2, i) = entity.z;
          } else if constexpr (std::is_same_v<T, Molecule>) {
            Vec3 O = {0, 0, 0};
            switch (o) {
            case Molecule::Cartesian:
              throw std::runtime_error(
                  "Can't use Cartesian origin for NeighbourList");
            case Molecule::Centroid:
              O = entity.centroid();
              break;
            case Molecule::CentreOfMass:
              O = entity.centre_of_mass();
              break;
            }
            cart_pos(0, i) = O.x();
            cart_pos(1, i) = O.y();
            cart_pos(2, i) = O.z();
          }
        },
        ent);
  }
  UnitCell uc = this->unit_cell();
  if (uc.dummy()) {
    Vec3 min_vals = cart_pos.rowwise().minCoeff();
    Mat3N shifted_cart_pos = cart_pos.colwise() - min_vals;
    Mat3N frac_pos = uc.to_fractional(shifted_cart_pos);
    nlp = NeighbourListPacket(entities, cart_pos, frac_pos);
  } else {
    Mat3N frac_pos = uc.to_fractional(cart_pos);
    frac_pos = frac_pos.array() - frac_pos.array().floor();
    Mat3N wrapped_cart_pos = uc.to_cartesian(frac_pos);
    nlp = NeighbourListPacket(entities, wrapped_cart_pos, frac_pos);
  }
  return nlp;
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

  Entities ents = this->get_entities();
  NeighbourListPacket nlp = this->get_neighpack_from_entities(ents);
  m_topo_neigh_list.update(nlp);

  // TODO: allow user input tolerance
  double tol = 0.4;
  BondGraph &bg = this->m_bond_graph;
  // TODO: make neighbourcallback func more generic
  //  to allow passing Atom instead of Entity to
  //  prevent the std::visit call every iteration
  NeighbourCallback func = [&ents, &bg, tol](const Entity &ent1,
                                             const Entity &ent2, double rsq) {
    const auto &ntt1 = ents[ent1.idx];
    const auto &ntt2 = ents[ent2.idx];
    std::visit(
        [rsq, &bg, tol](const auto &e1, const auto &e2) {
          using T1 = std::decay_t<decltype(e1)>;
          using T2 = std::decay_t<decltype(e2)>;
          if constexpr (std::is_same_v<T1, Atom> && std::is_same_v<T2, Atom>) {
            std::optional<Bond> bond = e1.is_bonded_with_rsq(e2, rsq, tol);
            if (bond) {
              bg.add_edge(e1.index, e2.index, *bond);
            }
          }
        },
        ntt1, ntt2);
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
