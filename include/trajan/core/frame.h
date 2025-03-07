#pragma once
#include <trajan/core/atom.h>
#include <trajan/core/unit_cell.h>
#include <vector>

namespace trajan::core {
class Frame {
public:
  void set_atoms(std::vector<Atom> &atoms);
  std::vector<Atom> get_atoms() const { return m_atoms; };
  Mat3N get_atom_positions() const { return m_atom_positions; };

  void set_uc(UnitCell &uc) { m_uc = uc; };
  UnitCell get_uc() { return m_uc; }

private:
  static constexpr int EMPTY = -1;
  std::vector<Atom> m_atoms;
  Mat3N m_atom_positions;
  UnitCell m_uc;
  int m_num_atoms = EMPTY;
  double m_time;
  int m_index = 0;
};
}; // namespace trajan::core
