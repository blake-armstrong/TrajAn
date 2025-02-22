#pragma once
#include <trajan/core/atom.h>
#include <trajan/core/unit_cell.h>
#include <vector>

const int EMPTY = -1;

namespace trajan::core {
class Frame {
public:
  void set_atoms(std::vector<Atom> &atoms);
  std::vector<Atom> get_atoms() const { return m_atoms; };

  void set_uc(UnitCell &uc) { m_uc = uc; };

private:
  std::vector<Atom> m_atoms;
  UnitCell m_uc;
  int m_num_atoms = EMPTY;
  double m_time;
  int m_index = 0;
};
}; // namespace trajan::core
