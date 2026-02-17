#include "trajan/core/topology.h"
#include <fmt/core.h>
#include <fstream>
#include <occ/crystal/unitcell.h>
#include <optional>
#include <stdexcept>
#include <trajan/core/atom.h>
// #include <trajan/core/element.h>
#include <trajan/core/frame.h>
#include <trajan/core/log.h>
// #include <trajan/core/unit_cell.h>
#include <trajan/core/units.h>
#include <trajan/io/pdb.h>

namespace trajan::io {

using occ::crystal::triclinic_cell;
using occ::crystal::UnitCell;
using trajan::core::Frame;
using trajan::units::radians;
using Atom = trajan::core::EnhancedAtom;

bool PDBHandler::_initialise() {
  switch (m_mode) {
  case Mode::Read:
    m_infile.open(this->file_path());
    return m_infile.is_open();
  case Mode::Write:
    m_outfile.open(this->file_path());
    return m_outfile.is_open();
  };
}

void PDBHandler::_finalise() {
  switch (m_mode) {
  case Mode::Read:
    if (m_infile.is_open()) {
      m_infile.close();
    }
    return;
  case Mode::Write:
    if (m_outfile.is_open()) {
      m_outfile.flush();
      m_outfile.close();
    }
    return;
  }
}

bool PDBHandler::parse_pdb(Frame &frame) {
  trajan::log::trace("Attempting to parse PDB file");
  std::vector<core::Atom> &atoms = frame.atoms();
  atoms.clear();
  bool atoms_set{false};
  core::AtomGraph atom_graph;
  std::string line;
  while (std::getline(m_infile, line)) {
    if (line.substr(0, 6) == "CRYST1") {
      double a, b, c, alpha, beta, gamma;
      char record_name[7], sg[12], z[5];
      int result =
          std::sscanf(line.c_str(), PDB_CRYST_FMT_READ.data(), record_name, &a,
                      &b, &c, &alpha, &beta, &gamma, sg, z);
      // sg[11] = '\0';
      // z[4] = '\0';
      trajan::log::trace("{}: a={:.4f} b={:.4f} c={:.4f} Angstroms",
                         record_name, a, b, c);
      trajan::log::trace("{}: α={:.4f} β={:.4f} γ={:.4f} degrees", record_name,
                         alpha, beta, gamma);
      UnitCell uc = triclinic_cell(a, b, c, radians(alpha), radians(beta),
                                   radians(gamma));
      frame.set_unit_cell(uc);
      continue;
    }
    if (line.substr(0, 4) == "ATOM" || line.substr(0, 6) == "HETATM") {
      Atom atom;
      char record_name[7], name_buffer[5], alt_loc[2], res_name[4], ins_code[2],
          tmp[12];
      char chain_id[2], element_buffer[3], charge[3];
      int serial, res_seq;
      float occupancy, temp_factor;
      int result = std::sscanf(line.c_str(), PDB_LINE_FMT_READ.data(),
                               record_name, &serial, tmp, name_buffer, alt_loc,
                               res_name, tmp, chain_id, &res_seq, ins_code, tmp,
                               &atom.x, &atom.y, &atom.z, &occupancy,
                               &temp_factor, tmp, element_buffer, charge);
      if (result < 16) {
        std::runtime_error(fmt::format("Failed to parse line: '{}'", line));
      }
      record_name[6] = '\0';
      name_buffer[4] = '\0';
      alt_loc[1] = '\0';
      res_name[3] = '\0';
      chain_id[1] = '\0';
      ins_code[1] = '\0';
      element_buffer[2] = '\0';
      charge[2] = '\0';
      atom.uindex = serial;
      atom.type = name_buffer;
      occ::util::trim(atom.type);
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
      atom.index = atoms.size();
      atom.molecule_type = res_name;
      occ::util::trim(atom.molecule_type);
      atom.umolecule_index = res_seq;
      atom.molecule_index = 0;
      atoms.push_back(atom);
      atoms_set = false;
      trajan::log::trace("{}: Creating atom: {}", record_name, atom.repr());
      atom_graph.add_vertex(trajan::core::AtomVertex{atom.index});
      continue;
    };
    if (line.substr(0, 6) == "CONECT") {
      if (!atoms_set) {
        frame.set_atoms(atoms);
        atoms_set = true;
      }
      char record_name[7];
      int ais[5];
      int result =
          std::sscanf(line.c_str(), PDB_CONECT_FMT_READ.data(), record_name,
                      &ais[0], &ais[1], &ais[2], &ais[3], &ais[4]);
      core::Mat3N cart_pos = frame.cart_pos();
      const auto &uc = frame.unit_cell();
      auto frac_pos =
          uc ? std::make_optional(uc->to_fractional(cart_pos)) : std::nullopt;
      size_t a_main = ais[0] - 1;
      for (int ai = 1; ai < result - 1; ai++) {
        core::Bond bond;
        size_t aj = ais[ai] - 1;
        bond.indices = {a_main, aj};
        if (uc) {
          if (!frac_pos) {
            throw std::runtime_error("No fractional coordinates");
          }
          if (a_main >= frac_pos->cols() || aj >= frac_pos->cols()) {
            throw std::runtime_error("Index out of bounds");
          }
          core::Vec3 frac_dist = frac_pos->col(a_main) - frac_pos->col(aj);
          frac_dist = frac_dist.array() - frac_dist.array().round();
          core::Vec3 cart_dist = uc->to_cartesian(frac_dist);
          bond.bond_length = std::sqrt(cart_dist.squaredNorm());
        } else {
          core::Vec3 cart_dist = cart_pos.col(a_main) - cart_pos.col(aj);
          bond.bond_length = std::sqrt(cart_dist.squaredNorm());
        }
        atom_graph.add_edge(a_main, aj, bond, true);
        trajan::log::trace(
            "CONECT: Creating bond between atom indices {:>5} {:>5}: "
            "{:>8.4f} Angstroms",
            a_main, aj, bond.bond_length);
        if (atoms[a_main].molecule_type != atoms[aj].molecule_type) {
          trajan::log::warn(
              "  atoms {:>5} {:>5} are bonded according to CONECT but not of "
              "the same residue type({} and {}) !",
              a_main, aj, atoms[a_main].molecule_type, atoms[aj].molecule_type);
        }
        if (atoms[a_main].molecule_type == atoms[aj].molecule_type) {
          if (atoms[a_main].molecule_index != atoms[aj].molecule_index) {
            trajan::log::warn(
                "  atoms {:>5} {:>5} are bonded according to CONECT but do not "
                "have the same residue index({} and {} of residue type{})!",
                a_main, aj, atoms[a_main].molecule_index,
                atoms[aj].molecule_index, atoms[a_main].molecule_type);
          }
        }
      }
      continue;
    }
    trajan::log::debug("Skipping line: '{}'", line);
  }
  if (!atoms_set) {
    frame.set_atoms(atoms);
  }
  // TODO: process atom graph and input residue names/indicies and check
  // agreement with the info from atom_graph.
  frame.set_atom_graph(atom_graph);
  if (atom_graph.num_edges() > 0) {
    frame.set_topology(core::Topology(atoms, atom_graph));
  }
  trajan::log::trace("Found {} atoms in PDB file", atoms.size());
  trajan::log::trace("Successfully parsed PDB file");
  return true;
}

bool PDBHandler::read_next_frame(Frame &frame) {
  if (m_has_read) {
    return false;
  }
  bool success = this->parse_pdb(frame);

  m_has_read = true;

  return success;
}

bool PDBHandler::write_next_frame(const Frame &frame) {
  const auto &atoms = frame.atoms();
  const auto &uc = frame.unit_cell();

  if (uc.has_value()) {
    trajan::log::trace("Writing unit cell info to PDB");
    const auto &unit_cell = uc.value();
    m_outfile << fmt::format("CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f} "
                             "P 1           1",
                             unit_cell.a(), unit_cell.b(), unit_cell.c(),
                             trajan::units::degrees(unit_cell.alpha()),
                             trajan::units::degrees(unit_cell.beta()),
                             trajan::units::degrees(unit_cell.gamma()))
              << std::endl;
  } else {
    trajan::log::trace("No unit cell info to write to PDB");
  }

  int i = 1;
  for (const auto &atom : atoms) {
    std::string line =
        fmt::format(PDB_LINE_FMT_WRITE.data(),
                    "ATOM",                 // field 1: 6 chars
                    i,                      // field 2: 5 digits
                    ' ',                    // field 3: 1 char
                    atom.type,              // field 4: 4 chars
                    ' ',                    // field 5: 1 char
                    "RES",                  // field 6: 3 chars (residue name)
                    ' ',                    // field 7: 1 char
                    'A',                    // field 8: 1 char (chain ID)
                    1,                      // field 9: 4 digits (resSeq)
                    ' ',                    // field 10: 1 char (iCode)
                    "",                     // field 11: 3 chars (altLoc?)
                    atom.x, atom.y, atom.z, // 3 coordinates, 8.3f
                    1.0,                    // occupancy
                    0.0,                    // tempFactor
                    "",                     // segment ID (10 chars)
                    atom.element.symbol(),  // 2 chars
                    ""                      // charge (2 chars)
        );
    trajan::log::trace("Writing PDB line: {}", line);
    m_outfile << line << std::endl;
    i++;
  }

  m_outfile << "ENDMDL" << std::endl;

  return true;
}

}; // namespace trajan::io
