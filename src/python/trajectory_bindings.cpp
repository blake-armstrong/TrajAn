#include "trajectory_bindings.h"
#include <nanobind/stl/filesystem.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <trajan/core/trajectory.h>

namespace nb = nanobind;
using namespace trajan::core;

void register_trajectory_bindings(nb::module_ &m) {
  nb::class_<Trajectory>(m, "Trajectory")
      .def(nb::init<>())
      .def("load_files", &Trajectory::load_files)
      .def("next_frame", &Trajectory::next_frame)
      .def("reset", &Trajectory::reset)
      .def_prop_ro("has_frames", &Trajectory::has_frames)
      .def_prop_ro("current_frame_index", &Trajectory::current_frame_index)
      .def_prop_ro("frame", &Trajectory::frame)
      .def_prop_ro("atoms", &Trajectory::atoms)
      .def_prop_ro("num_atoms", &Trajectory::num_atoms)
      .def_prop_ro("unit_cell", &Trajectory::unit_cell)
      .def("update_topology", nb::overload_cast<>(&Trajectory::update_topology))
      .def("update_topology", nb::overload_cast<std::vector<Atom> &>(
                                  &Trajectory::update_topology))
      .def_prop_ro("unit_cell_molecules", &Trajectory::unit_cell_molecules); //
}
