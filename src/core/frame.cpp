#include "trajan/core/topology.h"
#include "trajan/core/util.h"
#include <fmt/core.h>
#include <occ/crystal/unitcell.h>
#include <stdexcept>
#include <trajan/core/frame.h>
#include <trajan/core/log.h>

namespace trajan::core {

using occ::crystal::UnitCell;

void Frame::set_atoms(const std::vector<Atom> &atoms) {
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
}

void Frame::set_unit_cell(const UnitCell &unit_cell) {
  m_unit_cell = unit_cell;
}

const Mat3N Frame::cart_pos(std::optional<double> scale) const {
  if (m_num_atoms != m_atoms.size()) {
    trajan::log::warn("Number of atoms and atoms size in frame out of sync.");
  }
  Mat3N cart_pos(3, m_num_atoms);
  if (scale.has_value()) {
    const double scale_val = scale.value();
    for (size_t i = 0; i < m_num_atoms; i++) {
      const Atom &atom = m_atoms[i];
      cart_pos(0, i) = atom.x * scale_val;
      cart_pos(1, i) = atom.y * scale_val;
      cart_pos(2, i) = atom.z * scale_val;
    }
  } else {
    for (size_t i = 0; i < m_num_atoms; i++) {
      const Atom &atom = m_atoms[i];
      cart_pos(0, i) = atom.x;
      cart_pos(1, i) = atom.y;
      cart_pos(2, i) = atom.z;
    }
  }
  return cart_pos;
}

const occ::Vec Frame::cart_pos_flat(std::optional<double> scale) const {
  if (m_num_atoms != m_atoms.size()) {
    trajan::log::warn("Number of atoms and atoms size in frame out of sync.");
  }
  occ::Vec cart_pos(m_num_atoms * 3);
  if (scale.has_value()) {
    const double scale_val = scale.value();
    for (size_t i = 0; i < m_num_atoms; ++i) {
      Atom atom = m_atoms[i];
      cart_pos[3 * i] = atom.x * scale_val;
      cart_pos[3 * i + 1] = atom.y * scale_val;
      cart_pos[3 * i + 2] = atom.z * scale_val;
    }
  } else {
    for (size_t i = 0; i < m_num_atoms; ++i) {
      Atom atom = m_atoms[i];
      cart_pos[3 * i] = atom.x;
      cart_pos[3 * i + 1] = atom.y;
      cart_pos[3 * i + 2] = atom.z;
    }
  }
  return cart_pos;
}

std::vector<int> Frame::atomic_numbers() const {
  std::vector<int> atomic_numbers(m_atoms.size());
  for (int i = 0; i < m_atoms.size(); i++) {
    atomic_numbers[i] = m_atoms[i].atomic_number();
  }
  return atomic_numbers;
};

void Frame::update_atom_position(size_t idx, Vec3 &pos) {
  if (idx >= m_num_atoms) {
    throw std::runtime_error("Bad index.");
  }
  m_atoms[idx].set_position(pos);
}

void Frame::populate_angles(Topology &top) {
  const Mat3N cart_pos = this->cart_pos();
  Mat3N frac_pos;
  if (m_unit_cell) {
    frac_pos = m_unit_cell->to_fractional(cart_pos);
  }
  auto &angles = top.angles();
  trajan::log::debug("Populating angles theta");
  for (auto &angle : angles) {
    occ::Vec3 b1, b2;
    if (m_unit_cell) {
      b1 = frac_pos.col(angle.center_atom()) - frac_pos.col(angle.atom1());
      b1 = b1.array() - b1.array().round();
      b1 = m_unit_cell->to_cartesian(b1);
      b2 = frac_pos.col(angle.center_atom()) - frac_pos.col(angle.atom3());
      b2 = b2.array() - b2.array().round();
      b2 = m_unit_cell->to_cartesian(b2);
    } else {
      b1 = cart_pos.col(angle.center_atom()) - cart_pos.col(angle.atom1());
      b2 = cart_pos.col(angle.center_atom()) - cart_pos.col(angle.atom3());
    }
    angle.theta = trajan::util::angle_between(b1, b2);
  }
  auto &dihedrals = top.dihedrals();
  size_t d1, d2, d3, d4;
  trajan::log::debug("Populating dihedrals phi/distance");
  for (auto &dihedral : dihedrals) {
    occ::Vec3 b1, b2, b3;
    d1 = dihedral.atom_indices[0];
    d2 = dihedral.atom_indices[1];
    d3 = dihedral.atom_indices[2];
    d4 = dihedral.atom_indices[3];
    switch (dihedral.type) {
    case trajan::core::DihedralType::PROPER:
      b1 = frac_pos.col(d1) - frac_pos.col(d2);
      b1 = b1.array() - b1.array().round();
      b1 = m_unit_cell->to_cartesian(b1);
      b2 = frac_pos.col(d2) - frac_pos.col(d3);
      b2 = b2.array() - b2.array().round();
      b2 = m_unit_cell->to_cartesian(b2);
      b3 = frac_pos.col(d3) - frac_pos.col(d4);
      b3 = b3.array() - b3.array().round();
      b3 = m_unit_cell->to_cartesian(b3);
      dihedral.phi = trajan::util::dihedral_between(b1, b2, b3);
      break;
    case trajan::core::DihedralType::IMPROPER:
      b1 = frac_pos.col(d1) - frac_pos.col(d2);
      b1 = b1.array() - b1.array().round();
      b1 = m_unit_cell->to_cartesian(b1);
      b2 = frac_pos.col(d1) - frac_pos.col(d3);
      b2 = b2.array() - b2.array().round();
      b2 = m_unit_cell->to_cartesian(b2);
      b3 = frac_pos.col(d1) - frac_pos.col(d4);
      b3 = b3.array() - b3.array().round();
      b3 = m_unit_cell->to_cartesian(b3);
      dihedral.distance = trajan::util::out_of_plane_distance(b1, b2, b3);
      break;
    }
  }
}

} // namespace trajan::core
