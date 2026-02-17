#include "occ/core/units.h"
#include "trajan/core/frame.h"
#include "trajan/core/util.h"
#include <Eigen/Geometry>
#include <trajan/core/atomgraph.h>
#include <trajan/core/topology.h>

namespace trajan::core {

Topology::Topology(const std::vector<Atom> &atoms, const AtomGraph &atom_graph)
    : m_atom_graph(atom_graph), m_atoms(atoms) {

  using VD = AtomGraph::VertexDescriptor;
  using ED = AtomGraph::EdgeDescriptor;
  std::vector<std::tuple<VD, VD, ED>> result;
  ankerl::unordered_dense::set<VD> visited;

  auto collect_edges = [&](const VD &current, const VD &predecessor,
                           const ED &edge_desc) {
    visited.insert(current);
    if (current == predecessor)
      return;
    auto bond_pair = this->make_bond_pair(current, predecessor);
    const Bond &bond = atom_graph.edge(edge_desc);
    Bond forward_bond(bond.bond_length, bond_pair);
    m_bond_storage[bond_pair] = forward_bond;
  };

  for (const auto &[vertex_desc, vertex_data] : atom_graph.vertices()) {
    if (visited.contains(vertex_desc))
      continue;
    atom_graph.breadth_first_traversal_with_edge(vertex_desc, collect_edges);
  }

  this->generate_all_from_bonds();
}

Topology::Topology(const std::vector<Atom> &atoms) : m_atoms(atoms) {

  for (int i = 0; i < atoms.size(); i++) {
    m_atom_graph.add_vertex(trajan::core::AtomVertex{i});
  }
  for (int i = 0; i < atoms.size(); i++) {
    for (int j = i + 1; j < atoms.size(); j++) {
      auto bond_opt = atoms[i].is_bonded(atoms[j]);
      if (!bond_opt.has_value()) {
        continue;
      }
      Bond bond = bond_opt.value();
      auto bond_pair = this->make_bond_pair(i, j);
      bond.indices = bond_pair;

      m_atom_graph.add_edge(i, j, bond, true);
      m_bond_storage[bond_pair] = bond;
    }
  }
  this->generate_all_from_bonds();
}

void Topology::add_bond(size_t atom1, size_t atom2, double bond_length) {
  const auto &[p1, p2] = this->make_bond_pair(atom1, atom2);

  if (m_bond_storage.find({p1, p2}) != m_bond_storage.end()) {
    trajan::log::warn(
        fmt::format("Bond ({} {}) already in Topology.", atom1, atom2));
    return;
  }

  Bond bond(bond_length);
  bond.indices = {p1, p2};

  if (m_atom_graph.vertex(atom1) == m_atom_graph.vertices().end()) {
    m_atom_graph.add_vertex(AtomVertex{static_cast<int>(p1)});
  }
  if (m_atom_graph.vertex(atom2) == m_atom_graph.vertices().end()) {
    m_atom_graph.add_vertex(AtomVertex{static_cast<int>(p2)});
  }
  m_atom_graph.add_edge(p1, p2, bond, true);
  m_bond_storage[{p1, p2}] = bond;
}

// void Topology::remove_bond(size_t atom1, size_t atom2) {
//   auto bond_pair = this->make_bond_pair(atom1, atom2);
//
//   if (m_bond_storage.find(bond_pair) == m_bond_storage.end()) {
//     trajan::log::warn(
//         fmt::format("Bond ({} {}) not in Topology.", atom1, atom2));
//     return;
//   }
//
//   m_atom_graph.remove_edge(atom1, atom2);
//   m_bond_storage.erase(bond_pair);
// }

bool Topology::has_bond(size_t atom1, size_t atom2) const {
  auto bond_pair = this->make_bond_pair(atom1, atom2);
  return m_bond_storage.find(bond_pair) != m_bond_storage.end();
}

void Topology::clear_bonds() {
  // m_atom_graph.clear_edges();
  m_atom_graph = AtomGraph();
  m_bond_storage.clear();
  this->clear_angles();
  this->clear_dihedrals();
}

std::vector<Bond> Topology::get_bonds(bool bidrectional) const {
  std::vector<Bond> bonds;
  bonds.reserve(m_bond_storage.size());

  auto f = [&bonds](const auto &bond, const auto &bond_pair) {
    bonds.push_back(bond);
  };
  if (bidrectional) {
    auto f = [&bonds](const auto &bond, const auto &bond_pair) {
      bonds.push_back(bond);
      bonds.push_back(
          Bond(bond.bond_length, {bond_pair.second, bond_pair.first}));
    };
  }
  for (const auto &[bond_pair, bond] : m_bond_storage) {
    f(bond, bond_pair);
  }

  return bonds;
}

std::optional<Bond> Topology::get_bond(size_t atom1, size_t atom2) const {
  auto bond_pair = this->make_bond_pair(atom1, atom2);
  auto it = m_bond_storage.find(bond_pair);

  if (it != m_bond_storage.end()) {
    return it->second;
  }

  return std::nullopt;
}

bool Topology::update_bond(size_t atom1, size_t atom2,
                           const Bond &updated_bond) {
  auto bond_pair = this->make_bond_pair(atom1, atom2);
  auto it = m_bond_storage.find(bond_pair);

  if (it != m_bond_storage.end()) {
    Bond corrected_bond = updated_bond;
    corrected_bond.indices = bond_pair;
    it->second = corrected_bond;
    m_atom_graph.add_edge(atom1, atom2, corrected_bond, true);
    return true;
  }

  return false;
}

void Topology::add_angle(size_t atom1, size_t centre, size_t atom3) {
  Angle angle(atom1, centre, atom3);
  if (m_angle_set.find(angle) != m_angle_set.end()) {
    trajan::log::warn(fmt::format("Angle ({} {} {}) already in Topology.",
                                  atom1, centre, atom3));
    return;
  }
  m_angles.push_back(angle);
  this->update_angle_structures();
}

void Topology::remove_angle(size_t atom1, size_t centre, size_t atom3) {
  Angle target(atom1, centre, atom3);
  if (m_angle_set.find(target) == m_angle_set.end()) {
    trajan::log::warn(
        fmt::format("Angle ({} {} {}) not in Topology.", atom1, centre, atom3));
    return;
  }
  m_angles.erase(std::remove(m_angles.begin(), m_angles.end(), target),
                 m_angles.end());
  this->update_angle_structures();
}

bool Topology::has_angle(size_t atom1, size_t centre, size_t atom3) const {
  Angle target(atom1, centre, atom3);
  return m_angle_set.find(target) != m_angle_set.end();
}

void Topology::clear_angles() {
  m_angles.clear();
  m_angle_set.clear();
}

void Topology::add_dihedral(size_t atom1, size_t atom2, size_t atom3,
                            size_t atom4, DihedralType type) {
  Dihedral dihedral(atom1, atom2, atom3, atom4, type);
  if (m_dihedral_set.find(dihedral) != m_dihedral_set.end()) {
    trajan::log::warn(
        fmt::format("Dihedral ({} {} {} {}, {}) already in Topology.", atom1,
                    atom2, atom3, atom4, dihedral_type_to_string(type)));
    return;
  }
  m_dihedrals.push_back(dihedral);
  this->update_dihedral_structures();
}

void Topology::remove_dihedral(size_t atom1, size_t atom2, size_t atom3,
                               size_t atom4) {
  Dihedral target(atom1, atom2, atom3, atom4);
  if (m_dihedral_set.find(target) == m_dihedral_set.end()) {
    trajan::log::warn(fmt::format("Dihedral ({} {} {} {}) not in Topology.",
                                  atom1, atom2, atom3, atom4));
    return;
  }
  m_dihedrals.erase(std::remove(m_dihedrals.begin(), m_dihedrals.end(), target),
                    m_dihedrals.end());
  this->update_dihedral_structures();
}

bool Topology::has_dihedral(size_t atom1, size_t atom2, size_t atom3,
                            size_t atom4) const {
  Dihedral target1(atom1, atom2, atom3, atom4);
  bool test1 = m_dihedral_set.find(target1) != m_dihedral_set.end();
  Dihedral target2(atom4, atom3, atom2, atom1);
  bool test2 = m_dihedral_set.find(target2) != m_dihedral_set.end();
  return test1 || test2;
}

void Topology::clear_dihedrals() {
  m_dihedrals.clear();
  m_dihedral_set.clear();
}

void Topology::generate_angles_from_bonds() {
  this->clear_angles();

  const auto &adj_list = m_atom_graph.adjacency_list();
  for (const auto &[atom_idx, neighbours] : adj_list) {
    auto angles = this->find_angles_around_atom(atom_idx);
    for (const auto &angle : angles) {
      m_angles.push_back(angle);
    }
  }

  this->update_angle_structures();
}

void Topology::generate_proper_dihedrals_from_bonds() {
  m_dihedrals.erase(std::remove_if(m_dihedrals.begin(), m_dihedrals.end(),
                                   [](const Dihedral &d) {
                                     return d.type == DihedralType::PROPER;
                                   }),
                    m_dihedrals.end());

  std::vector<Dihedral> new_dihedrals;

  for (const auto &[bond_pair, bond] : m_bond_storage) {
    size_t atom1 = bond_pair.first;
    size_t atom2 = bond_pair.second;

    auto dihedrals = this->find_proper_dihedrals_for_bond(atom1, atom2);
    new_dihedrals.insert(new_dihedrals.end(), dihedrals.begin(),
                         dihedrals.end());
  }

  std::sort(new_dihedrals.begin(), new_dihedrals.end());
  new_dihedrals.erase(std::unique(new_dihedrals.begin(), new_dihedrals.end()),
                      new_dihedrals.end());

  for (const auto &dihedral : new_dihedrals) {
    m_dihedrals.push_back(dihedral);
  }

  this->update_dihedral_structures();
}

void Topology::generate_improper_dihedrals_from_bonds() {
  std::vector<Dihedral> improper_dihedrals;
  const auto &adj_list = m_atom_graph.adjacency_list();

  for (const auto &[atom_idx, neighbours] : adj_list) {
    if (neighbours.size() != 3) {
      continue;
    }
    auto impropers = this->find_improper_dihedrals_around_atom(atom_idx);
    improper_dihedrals.insert(improper_dihedrals.end(), impropers.begin(),
                              impropers.end());
  }

  for (const auto &dihedral : improper_dihedrals) {
    m_dihedrals.push_back(dihedral);
  }

  this->update_dihedral_structures();
}

void Topology::generate_cyclic_structures(int max_cycle_size) const {
  m_atom_graph.generate_cycles(max_cycle_size);
}

void Topology::generate_all_from_bonds() {
  this->generate_angles_from_bonds();
  this->generate_proper_dihedrals_from_bonds();
  this->generate_improper_dihedrals_from_bonds();
  this->generate_molecules();
}

std::vector<size_t> Topology::get_bonded_atoms(size_t atom_idx) const {
  const auto &neighbours = m_atom_graph.neighbors(atom_idx);
  if (neighbours == m_atom_graph.adjacency_list().end()) {
    return {};
  }
  std::vector<size_t> bonded;
  bonded.reserve(neighbours->second.size());
  for (const auto &[neighbour_id, bond] : neighbours->second) {
    bonded.push_back(neighbour_id);
  }
  return bonded;
}

std::vector<Bond> Topology::get_bonds_involving_atom(size_t atom_idx) const {
  std::vector<Bond> bonds;

  for (const auto &[bond_pair, bond] : m_bond_storage) {
    if (bond_pair.first == atom_idx || bond_pair.second == atom_idx) {
      bonds.push_back(bond);
    }
  }

  return bonds;
}

std::vector<size_t> Topology::get_atoms_at_distance(size_t atom_idx,
                                                    size_t distance) const {
  if (distance == 0)
    return {atom_idx};
  if (distance == 1)
    return this->get_bonded_atoms(atom_idx);

  ankerl::unordered_dense::set<size_t> current_level = {atom_idx};
  ankerl::unordered_dense::set<size_t> visited = {atom_idx};

  for (size_t d = 1; d < distance; ++d) {
    ankerl::unordered_dense::set<size_t> next_level;
    for (size_t atom : current_level) {
      auto neighbours = this->get_bonded_atoms(atom);
      for (size_t neighbour : neighbours) {
        if (!visited.contains(neighbour)) {
          next_level.insert(neighbour);
          visited.insert(neighbour);
        }
      }
    }
    current_level = std::move(next_level);
  }

  ankerl::unordered_dense::set<size_t> final_level;
  for (size_t atom : current_level) {
    auto neighbours = this->get_bonded_atoms(atom);
    for (size_t neighbour : neighbours) {
      if (!visited.contains(neighbour)) {
        final_level.insert(neighbour);
      }
    }
  }

  return std::vector<size_t>(final_level.begin(), final_level.end());
}

void Topology::generate_molecules() {

  auto components = m_atom_graph.connected_components();
  ankerl::unordered_dense::map<size_t, std::vector<AtomGraph::VertexDescriptor>>
      grouped;
  for (const auto &[vertex, component_id] : components) {
    grouped[component_id].push_back(vertex);
  }
  m_molecules.clear();
  m_atom_to_molecule.clear();
  using ankerl::unordered_dense::map;
  map<std::string, size_t> molecule_counts;
  size_t unique_molecule_counter{0};
  for (const auto &[component_id, vertices] : grouped) {
    std::vector<Atom> atoms;
    std::vector<std::string> residue_name;
    std::vector<size_t> residue_index;
    for (const auto vertex : vertices) {
      const auto &a = m_atoms[vertex];
      atoms.push_back(a);
      residue_name.push_back(a.molecule_type);
      residue_index.push_back(a.umolecule_index);
    }
    bool same_res_name = trajan::util::all_same(residue_name);
    if (!same_res_name) {
      trajan::log::warn("not all atoms have the same residue type in molecule "
                        "extracted from topology:\n {}",
                        trajan::util::format_vector(residue_name));
    }
    if (!trajan::util::all_same(residue_index) && same_res_name) {
      trajan::log::warn(
          "not all atoms have the same residue index in same "
          "residue name ({}) in molecule extracted from topology:\n {}",
          residue_name[0], trajan::util::format_vector(residue_name));
    }
    Molecule molecule(atoms);
    molecule.index = m_molecules.size();
    molecule.uindex = residue_index[0];
    molecule.utype = residue_name[0];
    molecule.type = fmt::format("M{}", unique_molecule_counter + 1);
    for (const auto &m : m_molecules) {
      if (molecule.is_comparable_to(m)) {
        molecule.type = m.type;
        break;
      }
      unique_molecule_counter++;
      molecule.type = fmt::format("M{}", unique_molecule_counter + 1);
    }
    molecule_counts[molecule.type]++;
    molecule.subindex = molecule_counts[molecule.type];
    m_molecules.emplace_back(atoms);
    for (const auto &a : atoms) {
      m_atom_to_molecule[a.index] = molecule.index;
    }
  }
}

size_t Topology::num_proper_dihedrals() const {
  return std::count_if(
      m_dihedrals.begin(), m_dihedrals.end(),
      [](const Dihedral &d) { return d.type == DihedralType::PROPER; });
}

size_t Topology::num_improper_dihedrals() const {
  return std::count_if(
      m_dihedrals.begin(), m_dihedrals.end(),
      [](const Dihedral &d) { return d.type == DihedralType::IMPROPER; });
}

void Topology::update_angle_structures() {
  m_angle_set.clear();
  for (const auto &angle : m_angles) {
    m_angle_set.insert(angle);
  }
}

void Topology::update_dihedral_structures() {
  m_dihedral_set.clear();
  for (const auto &dihedral : m_dihedrals) {
    m_dihedral_set.insert(dihedral);
  }
}

std::vector<Angle> Topology::find_angles_around_atom(size_t centre_atom) const {
  std::vector<Angle> angles;
  auto neighbours = this->get_bonded_atoms(centre_atom);

  if (neighbours.size() < 2) {
    return angles;
  }

  for (size_t i = 0; i < neighbours.size(); i++) {
    for (size_t j = i + 1; j < neighbours.size(); j++) {
      Angle angle(neighbours[i], centre_atom, neighbours[j]);
      angles.push_back(angle);
    }
  }

  return angles;
}

std::vector<Dihedral>
Topology::find_proper_dihedrals_for_bond(size_t atom1, size_t atom2) const {
  std::vector<Dihedral> dihedrals;

  auto neighbours1 = this->get_bonded_atoms(atom1);
  auto neighbours2 = this->get_bonded_atoms(atom2);

  neighbours1.erase(std::remove(neighbours1.begin(), neighbours1.end(), atom2),
                    neighbours1.end());
  neighbours2.erase(std::remove(neighbours2.begin(), neighbours2.end(), atom1),
                    neighbours2.end());

  for (size_t n1 : neighbours1) {
    for (size_t n2 : neighbours2) {
      Dihedral dihedral(n1, atom1, atom2, n2, DihedralType::PROPER);
      dihedrals.push_back(dihedral);
    }
  }

  return dihedrals;
}

std::vector<Dihedral>
Topology::find_improper_dihedrals_around_atom(size_t centre_atom) const {
  std::vector<Dihedral> impropers;
  auto neighbours = get_bonded_atoms(centre_atom);

  if (neighbours.size() == 3) {
    Dihedral improper(centre_atom, neighbours[0], neighbours[1], neighbours[2],
                      DihedralType::IMPROPER);
    impropers.push_back(improper);
  }

  return impropers;
}

bool Topology::validate_topology() const { return check_issues().empty(); }

std::vector<std::string> Topology::check_issues() const {
  std::vector<std::string> issues;

  for (const auto &angle : m_angles) {
    if (angle.atom1() == angle.center_atom() ||
        angle.atom3() == angle.center_atom() ||
        angle.atom1() == angle.atom3()) {
      issues.push_back(fmt::format("Invalid angle: {} {} {}", angle.atom1(),
                                   angle.center_atom(), angle.atom3()));
    }
  }

  for (const auto &dihedral : m_dihedrals) {
    ankerl::unordered_dense::set<size_t> unique_atoms(
        dihedral.atom_indices.begin(), dihedral.atom_indices.end());
    if (unique_atoms.size() != 4) {
      issues.push_back(
          fmt::format("Invalid dihedral: {} {} {} {}", dihedral.atom_indices[0],
                      dihedral.atom_indices[1], dihedral.atom_indices[2],
                      dihedral.atom_indices[3]));
    }
  }

  return issues;
}

void Topology::print_summary() const {
  size_t width = 20;
  trajan::log::info("Topology summary:");
  trajan::log::info("  {:>{}}: {:>10}", "atoms", width, this->num_atoms());
  trajan::log::info("  {:>{}}: {:>10}", "bonds", width, this->num_bonds());
  trajan::log::info("  {:>{}}: {:>10}", "angles", width, this->num_angles());
  trajan::log::info("  {:>{}}: {:>10}", "proper dihedrals", width,
                    this->num_proper_dihedrals());
  trajan::log::info("  {:>{}}: {:>10}", "improper dihedrals", width,
                    this->num_improper_dihedrals());
}

void Topology::print_detailed() const {
  trajan::log::info("Atoms:");
  trajan::log::info("  {:>8}  molecule-type(input) molecule-index(input)", "i");
  for (size_t i = 0; i < m_atoms.size(); i++) {
    const auto &atom = m_atoms[i];
    trajan::log::info("  {:>8d} {:<6s} {:>6}", i, atom.molecule_type,
                      atom.umolecule_index);
  }

  trajan::log::info("Bonds:");
  trajan::log::info("  {:>8} {:>8} length(Angstroms)", "i", "j");
  for (const auto &[bond_pair, bond] : m_bond_storage) {
    trajan::log::info("  {:>8d} {:>8d} {:>.4f} ", bond_pair.first,
                      bond_pair.second, bond.bond_length);
  }

  trajan::log::info("Angles:");
  trajan::log::info("  {:>8} {:>8} {:>8} theta(degrees)", "i", "j", "k");
  for (const auto &angle : m_angles) {
    trajan::log::info("  {:>8d} {:>8d} {:>8d} {:>.4f}", angle.atom1(),
                      angle.center_atom(), angle.atom3(),
                      occ::units::degrees(angle.theta));
  }

  trajan::log::info("Proper dihedrals:");
  trajan::log::info("  {:>8} {:>8} {:>8} {:>8} phi(degreees)", "i", "j", "k",
                    "l");
  for (const auto &dihedral : m_dihedrals) {
    switch (dihedral.type) {
    case DihedralType::PROPER:

      trajan::log::info("  {:>8d} {:>8d} {:>8d} {:>8d} {:>.4f} ",
                        dihedral.atom_indices[0], dihedral.atom_indices[1],
                        dihedral.atom_indices[2], dihedral.atom_indices[3],
                        occ::units::degrees(dihedral.phi));
      break;
    case DihedralType::IMPROPER:
      break;
    }
  }
  trajan::log::info("Improper dihedrals:");
  trajan::log::info(
      "  {:>8} {:>8} {:>8} {:>8} out-of-plane-distance(Angstroms)", "i", "j",
      "k", "l");
  for (const auto &dihedral : m_dihedrals) {
    switch (dihedral.type) {
    case DihedralType::PROPER:
      break;
    case DihedralType::IMPROPER:
      trajan::log::info("  {:>8d} {:>8d} {:>8d} {:>8d} {:>.4f} ",
                        dihedral.atom_indices[0], dihedral.atom_indices[1],
                        dihedral.atom_indices[2], dihedral.atom_indices[3],
                        dihedral.distance);
      break;
    }
  }

  trajan::log::info("Per-atom connectivity:");
  for (size_t i = 0; i < m_atoms.size(); i++) {
    auto bonded = get_bonded_atoms(i);

    size_t angle_count =
        std::count_if(m_angles.begin(), m_angles.end(), [i](const Angle &a) {
          return a.atom1() == i || a.center_atom() == i || a.atom3() == i;
        });
    size_t dihedral_count = std::count_if(
        m_dihedrals.begin(), m_dihedrals.end(), [i](const Dihedral &d) {
          return d.atom_indices[0] == i || d.atom_indices[1] == i ||
                 d.atom_indices[2] == i || d.atom_indices[3] == i;
        });

    std::string bonded_str;
    for (size_t j = 0; j < bonded.size(); j++) {
      bonded_str += std::to_string(bonded[j]);
      if (j + 1 < bonded.size())
        bonded_str += ", ";
    }

    trajan::log::info("  {:>8d}  bonds: {:<2d}  angles: {:<3d}  dihedrals: "
                      "{:<3d}  bonded to: [{}]",
                      i, bonded.size(), angle_count, dihedral_count,
                      bonded_str);
  }

  auto issues = check_issues();
  if (issues.empty()) {
    trajan::log::info("Validation: OK");
  } else {
    trajan::log::info("Validation: {} issue(s) found:", issues.size());
    for (const auto &issue : issues) {
      trajan::log::info("  [!] {}", issue);
    }
  }
}

std::string Topology::to_string() const {
  return fmt::format("Topology(bonds={}, angles={}, proper_dihedrals={}, "
                     "improper_dihedrals={})",
                     this->num_bonds(), this->num_angles(),
                     this->num_proper_dihedrals(),
                     this->num_improper_dihedrals());
}

} // namespace trajan::core
