#pragma once
#include <occ/core/atom.h>
#include <occ/core/element.h>
#include <occ/core/linear_algebra.h>
#include <occ/crystal/unitcell.h>
#include <optional>
#include <trajan/core/atomgraph.h>

namespace trajan::core {

using occ::Vec3;
using occ::core::Element;
using Bond = trajan::core::BondEdge;

static std::pair<Element, int>
convert_element_type(const std::string &element_type) {
  Element element = static_cast<Element>(element_type);
  return std::pair(element, element.atomic_number());
}

struct EnhancedAtom : public occ::core::Atom {
  Element element;
  std::string type, molecule_type = "UNK";
  int index, molecule_index;
  int uindex = 0, umolecule_index = 0;
  // size_t serial;

  EnhancedAtom() : Atom(0, 0.0, 0.0, 0.0), element(0) {};

  EnhancedAtom(int atomic_number, const Vec3 &pos)
      : Atom(atomic_number, pos[0], pos[1], pos[2]),
        element(static_cast<Element>(atomic_number)) {};

  EnhancedAtom(const std::string &element_type, const Vec3 &pos)
      : Atom(convert_element_type(element_type).second, pos[0], pos[1], pos[2]),
        element(convert_element_type(element_type).first) {};

  EnhancedAtom(int atomic_number, const Vec3 &pos, int index)
      : Atom(atomic_number, pos[0], pos[1], pos[2]),
        element(static_cast<Element>(atomic_number)), index(index) {};

  EnhancedAtom(const Vec3 &pos, const Element &element)
      : Atom(element.atomic_number(), pos[0], pos[1], pos[2]),
        element(element) {}

  EnhancedAtom(const Vec3 &pos, const Element &element, int index)
      : Atom(element.atomic_number(), pos[0], pos[1], pos[2]), element(element),
        index(index) {}

  EnhancedAtom(const Vec3 &pos, const std::string &element_str, int index)
      : Atom(convert_element_type(element_str).second, pos[0], pos[1], pos[2]),
        element(convert_element_type(element_str).second), index(index) {}

  inline std::optional<Bond> is_bonded(const EnhancedAtom &other,
                                       double bond_tolerance = 0.4) const {
    double rsq = this->square_distance(other);
    return this->is_bonded_with_sq_distance(other, rsq, bond_tolerance);
  }

  inline std::optional<Bond>
  is_bonded_with_rsq(const EnhancedAtom &other, double rsq,
                     double bond_tolerance = 0.4) const {
    return this->is_bonded_with_sq_distance(other, rsq, bond_tolerance);
  }

  inline const int atomic_number() const {
    return this->element.atomic_number();
  }

  inline bool operator==(const EnhancedAtom &rhs) const {
    return this->index == rhs.index;
  }

  inline const std::string repr() const {
    return fmt::format("Atom(index={}, type={}, element={}, x={:.3f}, "
                       "y={:.3f}, z={:.3f})",
                       index, type, element.symbol(), x, y, z);
  }

private:
  inline std::optional<Bond>
  is_bonded_with_sq_distance(const EnhancedAtom &other, double rsq,
                             double bond_tolerance) const {
    double distance = std::sqrt(rsq);
    double bond_threshold = element.covalent_radius() +
                            other.element.covalent_radius() + bond_tolerance;
    if (distance > bond_threshold) {
      return std::nullopt;
    }
    Bond bond(distance);
    return bond;
  }
};
}; // namespace trajan::core
