#include "core_bindings.h"
#include <Eigen/Core>
#include <nanobind/eigen/dense.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <trajan/core/atom.h>
#include <trajan/core/molecule.h>
#include <trajan/core/topology.h>

using namespace trajan::core;

nb::module_ register_core_bindings(nb::module_ &m) {
  using namespace nb::literals;

  m.doc() = "Trajan Core Module";

  // Bond structure
  nb::class_<Bond>(m, "Bond")
      .def(nb::init<>())
      .def(nb::init<double>(), "bond_length"_a)
      .def_rw("bond_length", &Bond::bond_length)
      .def_rw("idxs", &Bond::idxs)
      .def("__repr__", [](const Bond &b) {
        return fmt::format("Bond(length={:.3f}, idxs=({}, {}))", b.bond_length,
                           b.idxs.first, b.idxs.second);
      });

  // Angle structure
  nb::class_<Angle>(m, "Angle")
      .def(nb::init<>())
      .def(nb::init<size_t, size_t, size_t>(), "atom1"_a, "center"_a, "atom3"_a)
      .def_rw("atom_indices", &Angle::atom_indices)
      .def_rw("equilibrium_angle", &Angle::equilibrium_angle)
      .def_rw("force_constant", &Angle::force_constant)
      .def("center_atom", &Angle::center_atom)
      .def("atom1", &Angle::atom1)
      .def("atom3", &Angle::atom3)
      .def("__eq__", &Angle::operator==)
      .def("__repr__", [](const Angle &a) {
        return fmt::format("Angle({}, {}, {}, eq_angle={:.3f}, k={:.3f})",
                           a.atom1(), a.center_atom(), a.atom3(),
                           a.equilibrium_angle, a.force_constant);
      });

  // DihedralType enum
  nb::enum_<DihedralType>(m, "DihedralType")
      .value("PROPER", DihedralType::PROPER)
      .value("IMPROPER", DihedralType::IMPROPER)
      .export_values();

  // Dihedral structure
  nb::class_<Dihedral>(m, "Dihedral")
      .def(nb::init<>())
      .def(nb::init<size_t, size_t, size_t, size_t, DihedralType>(), "atom1"_a,
           "atom2"_a, "atom3"_a, "atom4"_a, "type"_a = DihedralType::PROPER)
      .def_rw("atom_indices", &Dihedral::atom_indices)
      .def_rw("type", &Dihedral::type)
      .def_rw("equilibrium_angle", &Dihedral::equilibrium_angle)
      .def_rw("force_constant", &Dihedral::force_constant)
      .def_rw("multiplicity", &Dihedral::multiplicity)
      .def("atom2", &Dihedral::atom2)
      .def("atom3", &Dihedral::atom3)
      .def("__eq__", &Dihedral::operator==)
      .def("__repr__", [](const Dihedral &d) {
        return fmt::format("Dihedral({}, {}, {}, {}, type={}, eq_angle={:.3f})",
                           d.atom_indices[0], d.atom_indices[1],
                           d.atom_indices[2], d.atom_indices[3],
                           dihedral_type_to_string(d.type),
                           d.equilibrium_angle);
      });

  // Topology class - the main class
  nb::class_<Topology>(m, "Topology")
      .def(nb::init<>())
      .def(nb::init<const BondGraph &>(), "bond_graph"_a)
      .def(nb::init<const std::vector<Atom> &>(), "atoms"_a)

      // Bond management
      .def("add_bond", &Topology::add_bond, "atom1"_a, "atom2"_a,
           "bond_length"_a = 0.0, "Add a bond between two atoms")
      .def("remove_bond", &Topology::remove_bond, "atom1"_a, "atom2"_a,
           "Remove a bond between two atoms")
      .def("has_bond", &Topology::has_bond, "atom1"_a, "atom2"_a,
           "Check if bond exists between two atoms")
      .def("clear_bonds", &Topology::clear_bonds, "Clear all bonds")
      .def("get_bonds", &Topology::get_bonds, "Get all bonds")
      .def("get_bond", &Topology::get_bond, "atom1"_a, "atom2"_a,
           "Get specific bond between two atoms")
      .def("update_bond", &Topology::update_bond, "atom1"_a, "atom2"_a,
           "bond"_a, "Update an existing bond")
      .def("get_bonds_involving_atom", &Topology::get_bonds_involving_atom,
           "atom_idx"_a, "Get all bonds involving a specific atom")

      // Angle management
      .def("add_angle", &Topology::add_angle, "atom1"_a, "center"_a, "atom3"_a,
           "Add an angle")
      .def("remove_angle", &Topology::remove_angle, "atom1"_a, "center"_a,
           "atom3"_a, "Remove an angle")
      .def("has_angle", &Topology::has_angle, "atom1"_a, "center"_a, "atom3"_a,
           "Check if angle exists")
      .def("clear_angles", &Topology::clear_angles, "Clear all angles")
      .def("get_angles", &Topology::get_angles,
           nb::rv_policy::reference_internal, "Get all angles")

      // Dihedral management
      .def("add_dihedral", &Topology::add_dihedral, "atom1"_a, "atom2"_a,
           "atom3"_a, "atom4"_a, "type"_a = DihedralType::PROPER,
           "Add a dihedral")
      .def("remove_dihedral", &Topology::remove_dihedral, "atom1"_a, "atom2"_a,
           "atom3"_a, "atom4"_a, "Remove a dihedral")
      .def("has_dihedral", &Topology::has_dihedral, "atom1"_a, "atom2"_a,
           "atom3"_a, "atom4"_a, "Check if dihedral exists")
      .def("clear_dihedrals", &Topology::clear_dihedrals, "Clear all dihedrals")
      .def("get_dihedrals", &Topology::get_dihedrals,
           nb::rv_policy::reference_internal, "Get all dihedrals")

      // Generation methods
      .def("generate_angles_from_bonds", &Topology::generate_angles_from_bonds,
           "Generate angles from existing bonds")
      .def("generate_proper_dihedrals_from_bonds",
           &Topology::generate_proper_dihedrals_from_bonds,
           "Generate proper dihedrals from existing bonds")
      .def("generate_improper_dihedrals_from_bonds",
           &Topology::generate_improper_dihedrals_from_bonds,
           "Generate improper dihedrals from existing bonds")
      .def("generate_all_from_bonds", &Topology::generate_all_from_bonds,
           "Generate all angles and dihedrals from existing bonds")

      // Graph queries
      .def("get_bonded_atoms", &Topology::get_bonded_atoms, "atom_idx"_a,
           "Get atoms bonded to a specific atom")
      .def("get_atoms_at_distance", &Topology::get_atoms_at_distance,
           "atom_idx"_a, "distance"_a, "Get atoms at specific graph distance")

      // Molecule extraction
      .def("extract_molecules", &Topology::extract_molecules,
           "Extract molecules as connected components")

      // Counts
      .def("num_bonds", &Topology::num_bonds, "Number of bonds")
      .def("num_angles", &Topology::num_angles, "Number of angles")
      .def("num_dihedrals", &Topology::num_dihedrals, "Number of dihedrals")
      .def("num_proper_dihedrals", &Topology::num_proper_dihedrals,
           "Number of proper dihedrals")
      .def("num_improper_dihedrals", &Topology::num_improper_dihedrals,
           "Number of improper dihedrals")

      // Validation
      .def("validate_topology", &Topology::validate_topology,
           "Validate topology consistency")
      .def("check_issues", &Topology::check_issues, "Check for topology issues")
      .def("print_summary", &Topology::print_summary, "Print topology summary")
      .def("to_string", &Topology::to_string,
           "Convert topology to string representation")

      // Access to underlying graph
      .def("get_bond_graph",
           static_cast<const BondGraph &(Topology::*)() const>(
               &Topology::get_bond_graph),
           nb::rv_policy::reference_internal,
           "Get underlying bond graph (const)")
      .def("get_bond_graph",
           static_cast<BondGraph &(Topology::*)()>(&Topology::get_bond_graph),
           nb::rv_policy::reference_internal,
           "Get underlying bond graph (mutable)")

      // String representation
      .def("__repr__", &Topology::to_string)
      .def("__str__", &Topology::to_string);

  // Utility functions
  m.def("dihedral_type_to_string", &dihedral_type_to_string, "type"_a,
        "Convert DihedralType enum to string");
  m.def("dihedral_type_from_string", &dihedral_type_from_string, "str"_a,
        "Convert string to DihedralType enum");
  return m;
}
