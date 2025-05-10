#pragma once
#include <trajan/core/atom.h>
#include <trajan/core/molecule.h>
#include <trajan/core/neigh.h>
#include <trajan/core/unit_cell.h>
#include <vector>

namespace trajan::core {

struct Frame {

public:
  enum UpdateFlags {
    UNINITIALISED = 0,
    UC_INIT = 1,
    ATOMS_INIT = 2,
    INITIALISED = 3
  };

  inline const UnitCell &unit_cell() const { return m_uc; }
  void set_uc(UnitCell &uc);

  inline const std::vector<Atom> &atoms() const { return m_atoms; }
  void set_atoms(const std::vector<Atom> &atoms);
  inline size_t num_atoms() const { return m_num_atoms; }

  inline const Mat3N &cart_pos() const { return m_cart_pos; }
  inline const Mat3N &wrapped_cart_pos() const { return m_wrapped_cart_pos; }
  inline const Mat3N &frac_pos() const { return m_frac_pos; }

  inline void set_frac_pos(Mat3N &frac_pos) { m_frac_pos = frac_pos; };
  inline void set_wrapped_cart_pos(Mat3N &wrapped_cart_pos) {
    m_wrapped_cart_pos = wrapped_cart_pos;
  };

private:
  UnitCell m_uc;
  std::vector<Atom> m_atoms;

  size_t m_num_atoms = 0;
  Mat3N m_cart_pos;
  Mat3N m_frac_pos;
  Mat3N m_wrapped_cart_pos;

  void update_positions();
  unsigned int m_positions_needs_update = UNINITIALISED;
};

}; // namespace trajan::core
