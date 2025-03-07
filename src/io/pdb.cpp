#include <fmt/core.h>
#include <fstream>
#include <stdexcept>
#include <trajan/core/atom.h>
#include <trajan/core/element.h>
#include <trajan/core/frame.h>
#include <trajan/core/log.h>
#include <trajan/io/file.h>

namespace trajan::io {

namespace core = trajan::core;

bool PDBHandler::read_frame(core::Frame &frame) {
  std::ifstream file(get_file_path());
  if (!file.is_open()) {
    throw std::runtime_error(
        fmt::format("Unable to open file for reading: '{}'", get_file_name()));
  }
  trajan::log::debug(
      fmt::format("Successfully opened file '{}'", get_file_name()));
  std::vector<core::Atom> atoms;
  std::string line;
  while (std::getline(file, line)) {
    if (line.substr(0, 6) == "CRYST1") {
      double a, b, c, alpha, beta, gamma;
      char record_name[7], sg[12], z[5];
      int result = std::sscanf(line.c_str(), PDB_CRYST_FMT.data(), record_name,
                               &a, &b, &c, &alpha, &beta, &gamma, sg, z);
      core::UnitCell uc(a, b, c, alpha, beta, gamma);
      frame.set_uc(uc);
      continue;
    }
    if (line.substr(0, 4) != "ATOM" && line.substr(0, 6) != "HETATM") {
      continue;
    };
    core::Atom atom;
    char record_name[7], name_buffer[5], altLoc, resName[4], i_code,
        element_buffer[3], charge[3], chainID;
    int serial, res_seq;
    double occupancy, temp_factor;

    // format string based on PDB standard
    int result = std::sscanf(line.c_str(), PDB_LINE_FMT.data(), record_name,
                             &serial, name_buffer, &altLoc, resName, &chainID,
                             &res_seq, &i_code, &atom.x, &atom.y, &atom.z,
                             &occupancy, &temp_factor, element_buffer, charge);

    if (result != 15) {
      std::runtime_error(fmt::format("Failed to parse line: '{}'", line));
    }

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
    atom.element = core::Element(element_identifier, exact);
    atoms.push_back(atom);
  }
  frame.set_atoms(atoms);
  return false;
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
