#include "trajan/io/file.h"
#include "trajan/core/trajectory.h"
#include <fstream>
#include <stdexcept>

namespace trajan::io {

namespace core = trajan::core;

core::Trajectory DCDHandler::read_trajectory(const std::string &filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Unable to open file for reading: " + filename);
  }

  core::Trajectory trajectory;
  // Read file contents and populate trajectory
  // ...

  return trajectory;
}

void DCDHandler::write_trajectory(const std::string &filename,
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
