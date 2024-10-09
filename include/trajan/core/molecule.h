#pragma once
#include <string>
#include <trajan/core/atom.h>
#include <vector>

namespace trajan::core {
class Molecule {
private:
  std::vector<Atom> atoms;
  std::string name;
  // Other properties...

public:
  // Methods to add atoms, get properties, etc.
};
}; // namespace trajan::core
