#include <fmt/core.h>
#include <fstream>
#include <stdexcept>
#include <trajan/core/atom.h>
#include <trajan/core/element.h>
#include <trajan/core/frame.h>
#include <trajan/core/log.h>
#include <trajan/core/unit_cell.h>
#include <trajan/core/units.h>
#include <trajan/io/pdb.h>

namespace trajan::io {

namespace core = trajan::core;
using trajan::units::radians;

bool PDBHandler::_initialise() {
  m_file.open(this->file_path());
  return m_file.is_open();
}

void PDBHandler::_finalise() {
  if (m_file.is_open()) {
    m_file.close();
  }
}

bool PDBHandler::parse_pdb(core::Frame &frame) {

  std::vector<core::Atom> atoms;
  std::string line;
  while (std::getline(m_file, line)) {
    if (line.substr(0, 6) == "CRYST1") {
      trajan::log::debug("Found unit cell from CRYST1 line in PDB.");
      double a, b, c, alpha, beta, gamma;
      char record_name[7], sg[12], z[5];
      int result =
          std::sscanf(line.c_str(), PDB_CRYST_FMT_READ.data(), record_name, &a,
                      &b, &c, &alpha, &beta, &gamma, sg, z);
      record_name[6] = '\0';
      sg[11] = '\0';
      z[4] = '\0';
      core::UnitCell uc = core::triclinic_cell(a, b, c, radians(alpha),
                                               radians(beta), radians(gamma));
      // core::UnitCell uc = trajan::core::cubic_cell(25);
      trajan::log::debug(fmt::format("uc: {}", uc.dummy()));
      frame.set_uc(uc);
      continue;
    }
    if (line.substr(0, 4) != "ATOM" && line.substr(0, 6) != "HETATM") {
      continue;
    };
    core::Atom atom;
    char record_name[7], name_buffer[5], alt_loc[2], res_name[4], ins_code[2],
        tmp[12];
    char chain_id[2], element_buffer[3], charge[3];
    int serial, res_seq;
    float occupancy, temp_factor;

    int result =
        std::sscanf(line.c_str(), PDB_LINE_FMT_READ.data(), record_name,
                    &serial, tmp, name_buffer, alt_loc, res_name, tmp, chain_id,
                    &res_seq, ins_code, tmp, &atom.x, &atom.y, &atom.z,
                    &occupancy, &temp_factor, tmp, element_buffer, charge);
    if (result != 15) {
      std::runtime_error(fmt::format("Failed to parse line: '{}'", line));
    }
    name_buffer[4] = '\0';
    alt_loc[1] = '\0';
    res_name[3] = '\0';
    chain_id[1] = '\0';
    ins_code[1] = '\0';
    element_buffer[2] = '\0';
    charge[2] = '\0';
    atom.serial = serial;
    atom.type = name_buffer;
    trajan::util::trim(atom.type);
    std::string element_identifier;
    bool exact;
    std::string element = element_buffer;
    if (!element.empty()) {
      exact = true;
      element_identifier = element;
    } else {
      exact = false;
      element_identifier = name_buffer;
    }
    atom.element = core::element::Element(element_identifier, exact);
    atom.index = atoms.size();
    atoms.push_back(atom);
  }
  frame.set_atoms(atoms);
  return true;
}

bool PDBHandler::read_next_frame(core::Frame &frame) {
  if (m_has_read) {
    return false;
  }
  bool success = this->parse_pdb(frame);

  m_has_read = true;

  return success;
}

// void PDBHandler::write_trajectory(const std::string &filename,
//                                   const core::Trajectory &trajectory) {
//   std::ofstream file(filename);
//   if (!file.is_open()) {
//     throw std::runtime_error("Unable to open file for writing: " + filename);
//   }

// Write trajectory to file
// For example:
// for (const auto &frame : trajectory.get_frames()) {
//   file << frame.getAtomCount() << "\n";
//   file << "Frame " << frame.getIndex() << "\n";
//   for (const auto &atom : frame.getAtoms()) {
//     file << std::setw(2) << atom.getElement() << " " << std::fixed
//          << std::setprecision(6) << atom.getX() << " " << atom.getY() << "
//          "
//          << atom.getZ() << "\n";
//   }
// }
}; // namespace trajan::io
