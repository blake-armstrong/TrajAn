#include <trajan/core/molecule.h>

namespace trajan::core {

Molecule::Molecule(const std::vector<Atom> &atoms)
    : m_atomic_numbers(atoms.size()), m_positions(3, atoms.size()),
      m_partial_charges(Vec::Zero(atoms.size())) {
  m_elements.reserve(atoms.size());
  for (size_t i = 0; i < atoms.size(); i++) {
    const auto &atom = atoms[i];
    m_elements.push_back(atom.element);
    m_atomic_numbers(i) = atom.element.atomic_number();
    // Internally store in angstroms
    m_positions(0, i) = atom.x;
    m_positions(1, i) = atom.y;
    m_positions(2, i) = atom.z;
  }
  m_name = chemical_formula(m_elements);
}

}; // namespace trajan::core
