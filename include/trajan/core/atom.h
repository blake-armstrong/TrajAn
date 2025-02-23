#pragma once
#include <trajan/core/element.h>
#include <trajan/core/linear_algebra.h>

namespace trajan::core {

struct Atom {
  // int atomic_number;
  double x, y, z;
  Element element;
  int index;
  Atom() : element(0), index(0) {}
  Atom(const Vec3 &pos, int atomic_number, int index)
      : x(pos[0]), y(pos[1]), z(pos[2]),
        element(static_cast<Element>(atomic_number)), index(index) {}
  Atom(const Atom &other, Vec3 shift)
      : x(other.x + shift.x()), y(other.y + shift.y()), z(other.z + shift.z()),
        element(other.element), index(other.index) {}

  inline Atom create_ghost(Vec3 shift) const { return Atom(*this, shift); }

  inline void update_position(const Vec3 &pos) {
    x = pos.x(), y = pos.y(), z = pos.z();
  }

  // convenience helper to convert this position into \a Vec3
  inline Vec3 position() const { return {x, y, z}; }

  inline int id() const { return index; }

  /// the square euclidean distance from another atom
  inline double square_distance(const Atom &other) const {
    double dx = other.x - x, dy = other.y - y, dz = other.z - z;
    return dx * dx + dy * dy + dz * dz;
  }

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
