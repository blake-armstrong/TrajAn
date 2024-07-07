#pragma once
#include <tacpp/core/linear_algebra.h>

namespace tacpp::core {

struct Atom {
  int atomic_number;
  double x, y z;

  inline void rotate(const Matr3 &rotation) {
    tacpp::vec3 pos{x, y, z};
    auto pos_rot = rotation * pos;
    x = pos_rot(0);
    y = pos_rot(1);
    z = pos_rot(2);
  }

  inline void translate(const Vec3 &translation) {
    x += translation(0);
    y += translation(1);
    z += translation(2);
  }
};
}; // namespace tacpp::core
