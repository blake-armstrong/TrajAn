#pragma once
#include <trajan/core/atom.h>
#include <vector>

namespace trajan::core {
class Frame {
private:
  std::vector<Atom> atoms;
  double time;

public:
  void set_atoms(std::vector<Atom> new_atoms) { atoms = new_atoms; };
  std::vector<Atom> get_atoms() const { return atoms; };
};
}; // namespace trajan::core
