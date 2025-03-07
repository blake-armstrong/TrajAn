#include <fmt/core.h>
#include <stdexcept>
// #include <trajan/core/atom.h>
#include <trajan/core/frame.h>
#include <trajan/core/log.h>

namespace trajan::core {
void Frame::set_atoms(std::vector<Atom> &atoms) {
  int num_atoms = atoms.size();
  if (num_atoms == 0) {
    throw std::runtime_error("No atoms!");
  }
  if (m_num_atoms == EMPTY) {
    m_num_atoms = num_atoms;
  }
  if (m_num_atoms != num_atoms) {
    trajan::log::debug(fmt::format(
        "Number of atoms has changed. Previously {} atoms, now {} atoms.",
        m_num_atoms, num_atoms));
  }
  m_num_atoms = num_atoms;
  m_atoms = atoms;
  Mat3N atoms_pos(3, num_atoms);
  for (int i = 0; i < m_num_atoms; ++i) {
    Atom atom = atoms[i];
    atoms_pos(0, i) = atom.x;
    atoms_pos(1, i) = atom.y;
    atoms_pos(2, i) = atom.z;
  }
  m_atom_positions = atoms_pos;
};

} // namespace trajan::core
