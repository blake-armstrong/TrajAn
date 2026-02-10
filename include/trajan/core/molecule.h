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
  std::string type;
  int index;
  std::optional<int> uindex;
  std::vector<Atom> enhanced_atoms;

  inline bool operator==(const EnhancedMolecule &rhs) const {
    return this->index == rhs.index;
  };

  EnhancedMolecule(const std::vector<Atom> &atoms);

  std::vector<int> atom_indices() const;

  inline const std::string repr() const {
    auto pos = this->center_of_mass();
    auto indices = this->atom_indices();
    auto elements = this->elements();
    std::vector<std::string> elements_str(elements.size());
    for (const auto &e : elements) {
      elements_str.push_back(e.symbol());
    }

    return fmt::format("Molecule(index={}, type={}, elements=[{}], atoms=[{}], "
                       "x={:.4f}, y={:.4f}, z={:.4f})",
                       index, type, trajan::util::format_vector(elements_str),
                       trajan::util::format_vector(indices), pos[0], pos[1],
                       pos[2]);
  }
};

}; // namespace trajan::core
