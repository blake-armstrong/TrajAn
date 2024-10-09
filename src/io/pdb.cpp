#include "trajan/core/trajectory.h"
#include "trajan/io/file.h"
#include <fmt/core.h>
#include <fstream>
#include <stdexcept>
#include <trajan/core/log.h>

namespace trajan::io {

namespace core = trajan::core;

core::Trajectory PDBHandler::read_trajectory(const std::string &filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Unable to open file for reading: " + filename);
  }

  std::vector<core::Atom> atoms;
  std::string line;
  while (std::getline(file, line)) {
    if (line.substr(0, 4) != "ATOM" && line.substr(0, 6) != "HETATM") {
      continue;
    };
    core::Atom atom;
    char record_name[7];

    // Format string based on PDB standard
    int result = std::sscanf(
        line.c_str(),
        "%6s%5d %4s%c%3s %c%4d%c   %8lf%8lf%8lf%6lf%6lf          %2s%2s",
        record_name, &atom.serial, atom.name, &atom.altLoc, atom.resName,
        &atom.chainID, &atom.res_seq, &atom.i_code, &atom.x, &atom.y, &atom.z,
        &atom.occupancy, &atom.temp_factor, atom.element, atom.charge);

    if (result == 15) {
      atoms.push_back(atom);
    } else {
      std::runtime_error(
          fmt::format("Warning: Failed to parse line: '{}'", line));
    }
  }
  core::Trajectory trajectory;
  core::Frame frame;
  frame.set_atoms(atoms);

  return trajectory;
}

void PDBHandler::write_trajectory(const std::string &filename,
                                  const core::Trajectory &trajectory) {
  std::ofstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Unable to open file for writing: " + filename);
  }

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
}
}; // namespace trajan::io
