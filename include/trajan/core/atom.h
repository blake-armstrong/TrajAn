#pragma once
#include <trajan/core/linear_algebra.h>

namespace trajan::core {

struct Atom {
  int atomic_number;
  double x, y, z;

  int serial;
  char name[5];
  char altLoc;
  char resName[4];
  char chainID;
  int res_seq;
  char i_code;
  double occupancy;
  double temp_factor;
  char element[3];
  char charge[3];

  inline void rotate(const trajan::Mat3 &rotation) {
    trajan::Vec3 pos{x, y, z};
    auto pos_rot = rotation * pos;
    x = pos_rot(0);
    y = pos_rot(1);
    z = pos_rot(2);
  }

  inline void translate(const trajan::Vec3 &translation) {
    x += translation(0);
    y += translation(1);
    z += translation(2);
  }
};
}; // namespace trajan::core
