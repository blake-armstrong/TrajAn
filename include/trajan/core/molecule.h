#pragma once
#include <string>
#include <trajan/core/atom.h>
#include <trajan/core/element.h>
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
}; // namespace trajan::core
