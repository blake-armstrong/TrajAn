#pragma once
#include <trajan/core/frame.h>
#include <trajan/core/molecule.h>
#include <vector>

namespace trajan::core {
class Trajectory {
private:
  std::vector<Frame> frames;
  std::vector<Molecule> molecules;
  // Other trajectory-wide properties...

public:
  // Methods to load frames, analyze data, etc.
};
}; // namespace trajan::core
