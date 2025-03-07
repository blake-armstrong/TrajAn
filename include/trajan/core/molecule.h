#pragma once
#include <string>
#include <trajan/core/atom.h>
#include <trajan/core/element.h>
#include <trajan/core/graph.h>
#include <trajan/core/linear_algebra.h>
#include <vector>

namespace trajan::core {

class Molecule {
public:
  inline explicit Molecule() {};

  /**
   * Construct a Molecule from a vector of Atom objects
   *
   * \param atoms std::vector<Atom> of length N atoms.
   *
   **/
  Molecule(const std::vector<Atom> &atoms);
  /**
   * The number of atoms in this molecule.
   *
   * \returns size size_t representing the number of atoms in this Molecule.
   *
   * Calculated based on the internal IVec of atomic numbers.
   */
  inline size_t size() const { return m_atomic_numbers.size(); }

private:
  std::vector<Atom> m_atoms;
  std::string m_name{""};
  std::vector<std::pair<size_t, size_t>> m_bonds;
  std::vector<Element> m_elements;
  IVec m_atomic_numbers;
  Mat3N m_positions;
  Vec m_partial_charges;
};

class MoleculeGraph : public Graph<Atom, Bond> {
public:
  MoleculeGraph(const std::vector<Atom> &atoms, double bond_tolerance = 0.4)
      : Graph<Atom, Bond>(atoms,
                          [bond_tolerance](const Atom &a1, const Atom &a2) {
                            return a1.is_bonded(a2, bond_tolerance);
                          }) {};
  NodeId get_node_id_from_node(const Atom &atom) const override {
    return atom.id();
  };
};

std::vector<Molecule> identify_molecules(const std::vector<Atom> &atoms,
                                         double bond_tolerance = 0.4);

}; // namespace trajan::core
