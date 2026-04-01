#pragma once
#include "trajan/core/util.h"
#include <occ/core/bondgraph.h>
#include <occ/core/molecule.h>
#include <string>
#include <trajan/core/atom.h>
#include <vector>

namespace trajan::core {

using Atom = trajan::core::EnhancedAtom;

class EnhancedMolecule : public occ::core::Molecule {
public:
  std::string type{"UNK"}, utype{"UNK"};
  int index, uindex{0}, subindex;
  std::vector<Atom> enhanced_atoms;

  inline bool operator==(const EnhancedMolecule &rhs) const {
    return this->index == rhs.index;
  };

  inline std::vector<std::string> element_symbols() const {
    std::vector<std::string> elements_str;
    elements_str.reserve(this->elements().size());
    for (const auto &e : this->elements()) {
      elements_str.push_back(e.symbol());
    }
    return elements_str;
  }

  // Sync the occ::Molecule internal position matrix from the current
  // enhanced_atom x/y/z coordinates (which are kept up to date per-frame).
  // occ::Molecule::positions() is read-only, so we use const_cast.
  // m_positions is stored in Angstroms (same units as enhanced_atoms).
  inline void sync_occ_positions() {
    auto &pos = const_cast<occ::Mat3N &>(this->positions());
    for (size_t i = 0; i < enhanced_atoms.size(); ++i) {
      pos(0, i) = enhanced_atoms[i].x;
      pos(1, i) = enhanced_atoms[i].y;
      pos(2, i) = enhanced_atoms[i].z;
    }
  }

  inline std::vector<std::string> atom_types() const {
    std::vector<std::string> atom_types;
    atom_types.reserve(enhanced_atoms.size());
    for (const auto &a : enhanced_atoms) {
      atom_types.push_back(a.type);
    }
    return atom_types;
  }

  EnhancedMolecule(const std::vector<Atom> &atoms);

  std::vector<int> atom_indices() const;
  inline const std::vector<Atom> &atoms() const { return enhanced_atoms; }
  inline std::vector<Atom> &atoms() { return enhanced_atoms; }

  inline const std::string repr() const {
    auto pos = this->center_of_mass();
    auto indices = this->atom_indices();
    auto elements = this->elements();
    std::vector<std::string> elements_str = this->element_symbols();

    return fmt::format("Molecule(index={}, type={}, elements=[{}], atoms=[{}], "
                       "x={:.4f}, y={:.4f}, z={:.4f})",
                       index, type, trajan::util::format_vector(elements_str),
                       trajan::util::format_vector(indices), pos[0], pos[1],
                       pos[2]);
  }
};

}; // namespace trajan::core
