#include "trajan/core/util.h"
#include <cmath>
#include <memory>
#include <occ/core/linear_algebra.h>
#include <occ/crystal/unitcell.h>
#include <trajan/io/xtc.h>

namespace trajan::io {

using occ::crystal::UnitCell;
using trajan::core::Frame;

bool XTCHandler::_initialise() {
  if (m_mode == Mode::Read) {
    m_infile.open(this->file_path(), std::ios::binary);
    if (!m_infile.is_open()) {
      return false;
    }
    m_infile.close();
    m_xtcreader = std::make_unique<XTCReader>(this->file_path());
    return true;
  } else {
    try {
      m_xtcwriter = std::make_unique<XTCWriter>(this->file_path());
    } catch (const std::exception &) {
      return false;
    }
    return m_xtcwriter->is_open();
  }
}

void XTCHandler::_finalise() {}

bool XTCHandler::read_next_frame(Frame &frame) {
  if (m_xtcreader->eot()) {
    return false;
  }
  m_xtcreader->next_frame();

  occ::Mat3 m;
  for (int row = 0; row < 3; ++row) {
    for (int col = 0; col < 3; ++col) {
      // convert to column major and from nm to ang
      m(col, row) = static_cast<double>(m_xtcreader->box[row][col]) * 10.0;
    }
  }
  {
    // Columns of m are the lattice vectors
    occ::Vec3 a_vec = m.col(0), b_vec = m.col(1), c_vec = m.col(2);
    double a = a_vec.norm(), b = b_vec.norm(), c = c_vec.norm();
    double cos_alpha = std::clamp(b_vec.dot(c_vec) / (b * c), -1.0, 1.0);
    double cos_beta  = std::clamp(a_vec.dot(c_vec) / (a * c), -1.0, 1.0);
    double cos_gamma = std::clamp(a_vec.dot(b_vec) / (a * b), -1.0, 1.0);
    double alpha = std::acos(cos_alpha), beta = std::acos(cos_beta),
           gamma = std::acos(cos_gamma);
    if (trajan::util::unitcell_is_reasonable(a, b, c, alpha, beta, gamma)) {
      UnitCell uc = occ::crystal::triclinic_cell(a, b, c, alpha, beta, gamma);
      frame.set_unit_cell(uc);
    }
  }

  size_t natoms = m_xtcreader->X.size() / 3;
  for (size_t i = 0; i < natoms; i++) {
    int iatom = i * 3;
    occ::Vec3 pos = {m_xtcreader->X[iatom] * 10.0,
                     m_xtcreader->X[iatom + 1] * 10.0,
                     m_xtcreader->X[iatom + 2] * 10.0};
    frame.update_atom_position(i, pos);
  }

  return true;
}

bool XTCHandler::write_next_frame(const Frame &frame) {
  size_t natoms = frame.num_atoms();
  uint32_t step = 0;
  float time = 0.0;
  occ::Mat3 m = occ::Mat3::Zero();
  if (frame.has_unit_cell()) {
    m = frame.unit_cell().value().direct();
  };
  std::array<std::array<float, 3>, 3> box;
  for (int col = 0; col < 3; ++col) {
    for (int row = 0; row < 3; ++row) {
      box[row][col] = m(col, row) / 10.0;
    }
  }
  occ::Mat3N mat = frame.cart_pos();
  occ::FMat3N fmat = (mat / 10.0).cast<float>();
  std::vector<float> coords(fmat.data(), fmat.data() + fmat.size());

  return m_xtcwriter->write_frame(natoms, step, time, box, coords.data());
}

} // namespace trajan::io
