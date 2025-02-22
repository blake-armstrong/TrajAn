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
};

} // namespace trajan::core
