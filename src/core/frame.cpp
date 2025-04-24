#include <fmt/core.h>
#include <trajan/core/frame.h>
#include <trajan/core/log.h>

namespace trajan::core {

void Frame::set_atoms(std::vector<Atom> &atoms) {
  int num_atoms = atoms.size();
  if (num_atoms == 0) {
    throw std::runtime_error("No atoms!");
  }
  if (m_num_atoms != num_atoms && m_num_atoms != 0) {
    trajan::log::debug(fmt::format(
        "Number of atoms has changed. Previously {} atoms, now {} atoms.",
        m_num_atoms, num_atoms));
  }
  m_num_atoms = num_atoms;
  m_atoms = atoms;
  Mat3N cart_pos(3, m_num_atoms);
  for (size_t i = 0; i < m_num_atoms; ++i) {
    Atom atom = m_atoms[i];
    cart_pos(0, i) = atom.x;
    cart_pos(1, i) = atom.y;
    cart_pos(2, i) = atom.z;
  }
  m_cart_pos = cart_pos;
  m_positions_needs_update |= ATOMS_INIT;
}
} // namespace trajan::core
