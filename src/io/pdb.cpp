#include "trajan/core/topology.h"
#include <fmt/core.h>
#include <fstream>
#include <occ/crystal/unitcell.h>
#include <optional>
#include <stdexcept>
#include <trajan/core/atom.h>
#include <trajan/core/frame.h>
#include <trajan/core/log.h>
#include <trajan/core/units.h>
#include <trajan/io/pdb.h>
#include <ankerl/unordered_dense.h>

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
    if (m_infile.is_open())
      m_infile.close();
    return;
  case Mode::Write:
    if (m_outfile.is_open()) {
      // Write CONECT records (topology, written once after all frames).
      for (const auto &cl : m_conect_lines)
        m_outfile << cl << std::endl;
      if (m_write_frame_count > 0)
        m_outfile << "END" << std::endl;
      m_outfile.flush();
      m_outfile.close();
    }
    return;
  }
}

bool PDBHandler::parse_pdb(Frame &frame) {
  trajan::log::trace("Attempting to parse PDB frame");
  std::vector<core::Atom> &atoms = frame.atoms();
  atoms.clear();

  // Carry forward unit cell from a previous frame if this frame has no CRYST1.
  if (m_cached_uc.has_value())
    frame.set_unit_cell(m_cached_uc.value());

  bool atoms_set{false};
  // Set true when ENDMDL is seen. We keep reading after ENDMDL to pick up
  // CONECT records that appear between ENDMDL and the next MODEL/ATOM/END.
  bool past_atoms{false};
  core::AtomGraph atom_graph;
  ankerl::unordered_dense::map<int, size_t> serial_to_index;
  core::Mat3N conect_cart_pos;
  std::optional<UnitCell> conect_uc;
  std::optional<core::Mat3N> conect_frac_pos;

  // Process a single line. Returns false when the frame boundary is reached.
  auto process = [&](const std::string &line) -> bool {
    const bool is_endmdl = line.size() >= 6 && line.substr(0, 6) == "ENDMDL";
    const bool is_end    = line.size() >= 3 && line.substr(0, 3) == "END" && !is_endmdl;
    const bool is_model  = line.size() >= 5 && line.substr(0, 5) == "MODEL";
    const bool is_atom   = line.size() >= 4 && (line.substr(0, 4) == "ATOM" ||
                                                 line.substr(0, 6) == "HETATM");

    if (is_end)
      return false; // hard stop

    if (is_endmdl) {
      past_atoms = true; // coordinate block done; keep reading for CONECT
      return true;
    }

    if (is_model) {
      if (!atoms.empty() || past_atoms) {
        // Next frame is starting — stash this line for the next call.
        m_pending_line = line;
        m_has_pending = true;
        return false;
      }
      return true; // MODEL before first ATOM — skip
    }

    if (is_atom) {
      if (past_atoms) {
        // New ATOM record after ENDMDL with no MODEL separator — next frame.
        m_pending_line = line;
        m_has_pending = true;
        return false;
      }
      Atom atom;
      char record_name[7], name_buffer[5], alt_loc[2], res_name[4],
           ins_code[2], tmp[12], chain_id[2], element_buffer[3], charge[3];
      int serial, res_seq;
      float occupancy, temp_factor;
      int result = std::sscanf(line.c_str(), PDB_LINE_FMT_READ.data(),
                               record_name, &serial, tmp, name_buffer, alt_loc,
                               res_name, tmp, chain_id, &res_seq, ins_code, tmp,
                               &atom.x, &atom.y, &atom.z, &occupancy,
                               &temp_factor, tmp, element_buffer, charge);
      if (result < 16)
        std::runtime_error(fmt::format("Failed to parse line: '{}'", line));
      record_name[6] = '\0';
      name_buffer[4] = '\0';
      alt_loc[1]     = '\0';
      res_name[3]    = '\0';
      chain_id[1]    = '\0';
      ins_code[1]    = '\0';
      element_buffer[2] = '\0';
      charge[2]      = '\0';
      atom.uindex = serial;
      atom.type = name_buffer;
      occ::util::trim(atom.type);
      std::string element = element_buffer;
      if (!element.empty()) {
        atom.element = core::Element(element, true);
      } else {
        atom.element = core::Element(std::string(name_buffer), false);
      }
      atom.index = atoms.size();
      atom.molecule_type = res_name;
      occ::util::trim(atom.molecule_type);
      atom.umolecule_index = res_seq;
      atom.molecule_index  = 0;
      serial_to_index[serial] = atoms.size();
      atoms.push_back(atom);
      atoms_set = false;
      trajan::log::trace("{}: Creating atom: {}", record_name, atom.repr());
      atom_graph.add_vertex(trajan::core::AtomVertex{atom.index});
      return true;
    }

    if (line.size() >= 6 && line.substr(0, 6) == "CRYST1") {
      double a, b, c, alpha, beta, gamma;
      char record_name[7], sg[12], z[5];
      std::sscanf(line.c_str(), PDB_CRYST_FMT_READ.data(), record_name,
                  &a, &b, &c, &alpha, &beta, &gamma, sg, z);
      trajan::log::trace("{}: a={:.4f} b={:.4f} c={:.4f} Angstroms",
                         record_name, a, b, c);
      UnitCell uc = triclinic_cell(a, b, c, radians(alpha), radians(beta),
                                   radians(gamma));
      m_cached_uc = uc;
      frame.set_unit_cell(uc);
      return true;
    }

    if (line.size() >= 6 && line.substr(0, 6) == "CONECT") {
      if (!atoms_set) {
        frame.set_atoms(atoms);
        atoms_set = true;
        conect_cart_pos = frame.cart_pos();
        conect_uc = frame.unit_cell();
        if (conect_uc)
          conect_frac_pos = conect_uc->to_fractional(conect_cart_pos);
      }
      auto read_col5 = [&](int col_start) -> int {
        if (col_start >= static_cast<int>(line.size()))
          return 0;
        int len = std::min(5, static_cast<int>(line.size()) - col_start);
        try { return std::stoi(line.substr(col_start, len)); }
        catch (...) { return 0; }
      };
      int main_serial = read_col5(6);
      if (main_serial <= 0)
        return true;
      auto main_it = serial_to_index.find(main_serial);
      if (main_it == serial_to_index.end()) {
        trajan::log::warn("CONECT: main serial {} not found, skipping", main_serial);
        return true;
      }
      size_t a_main = main_it->second;
      for (int col : {11, 16, 21, 26}) {
        int bonded_serial = read_col5(col);
        if (bonded_serial <= 0)
          continue;
        auto bonded_it = serial_to_index.find(bonded_serial);
        if (bonded_it == serial_to_index.end()) {
          trajan::log::warn("CONECT: bonded serial {} not found, skipping", bonded_serial);
          continue;
        }
        size_t aj = bonded_it->second;
        core::Bond bond;
        bond.indices = {static_cast<int>(a_main), static_cast<int>(aj)};
        if (conect_uc) {
          core::Vec3 frac_dist = conect_frac_pos->col(a_main) - conect_frac_pos->col(aj);
          frac_dist = frac_dist.array() - frac_dist.array().round();
          core::Vec3 cart_dist = conect_uc->to_cartesian(frac_dist);
          bond.bond_length = std::sqrt(cart_dist.squaredNorm());
        } else {
          core::Vec3 cart_dist = conect_cart_pos.col(a_main) - conect_cart_pos.col(aj);
          bond.bond_length = std::sqrt(cart_dist.squaredNorm());
        }
        atom_graph.add_edge(a_main, aj, bond, true);
        trajan::log::trace("CONECT: bond {:>5} {:>5}: {:>8.4f} A",
                           a_main, aj, bond.bond_length);
        if (atoms[a_main].molecule_type != atoms[aj].molecule_type)
          trajan::log::warn("  atoms {:>5} {:>5} bonded but different residue types ({} and {})",
                            a_main, aj, atoms[a_main].molecule_type, atoms[aj].molecule_type);
      }
      return true;
    }

    trajan::log::debug("Skipping line: '{}'", line);
    return true;
  };

  // Process the pending line from the previous call first.
  if (m_has_pending) {
    m_has_pending = false;
    if (!process(m_pending_line))
      goto done;
  }

  {
    std::string line;
    while (std::getline(m_infile, line)) {
      if (!process(line))
        break;
    }
  }

done:
  trajan::log::trace("Found {} atoms in PDB frame", atoms.size());
  if (atoms.empty()) {
    trajan::log::trace("No atoms found — end of PDB file");
    return false;
  }
  if (!atoms_set)
    frame.set_atoms(atoms);
  frame.set_atom_graph(atom_graph);
  if (atom_graph.num_edges() > 0)
    frame.set_topology(core::Topology(atoms, atom_graph));
  return true;
}

bool PDBHandler::read_next_frame(Frame &frame) {
  return this->parse_pdb(frame);
}

bool PDBHandler::write_next_frame(const Frame &frame) {
  const auto &atoms = frame.atoms();
  const auto &uc = frame.unit_cell();

  ++m_write_frame_count;

  // CRYST1 once at the top of the file.
  if (m_write_frame_count == 1) {
    if (uc.has_value()) {
      const auto &unit_cell = uc.value();
      m_outfile << fmt::format(
                       "CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f} "
                       "P 1           1",
                       unit_cell.a(), unit_cell.b(), unit_cell.c(),
                       trajan::units::degrees(unit_cell.alpha()),
                       trajan::units::degrees(unit_cell.beta()),
                       trajan::units::degrees(unit_cell.gamma()))
                << std::endl;
    }

    // Cache CONECT lines from the first frame that has topology.
    // They are written once at the end of the file (VMD-compatible layout).
    const auto &graph = frame.get_atom_graph();
    if (graph.num_edges() > 0) {
      const auto &adj = graph.adjacency_list();
      for (const auto &[vertex, neighbors] : adj) {
        if (neighbors.empty())
          continue;
        int main_serial = m_original_ids
                              ? atoms[vertex].uindex
                              : static_cast<int>(vertex) + 1;
        std::vector<int> bonded;
        bonded.reserve(neighbors.size());
        for (const auto &[nb, _] : neighbors)
          bonded.push_back(m_original_ids ? atoms[nb].uindex
                                          : static_cast<int>(nb) + 1);
        for (size_t i = 0; i < bonded.size(); i += 4) {
          std::string cl = fmt::format("CONECT{:5d}", main_serial);
          for (size_t j = i; j < std::min(i + 4, bonded.size()); ++j)
            cl += fmt::format("{:5d}", bonded[j]);
          m_conect_lines.push_back(cl);
        }
      }
    }
  }

  m_outfile << fmt::format("MODEL {:>8}", m_write_frame_count) << std::endl;

  int i = 1;
  for (const auto &atom : atoms) {
    int serial = m_original_ids ? atom.uindex : i;
    m_outfile << fmt::format(PDB_LINE_FMT_WRITE.data(),
                             "ATOM", serial, ' ', atom.type, ' ',
                             atom.molecule_type, ' ', 'A',
                             atom.umolecule_index, ' ', "",
                             atom.x, atom.y, atom.z,
                             1.0, 0.0, "",
                             atom.element.symbol(), "")
              << std::endl;
    i++;
  }

  m_outfile << "ENDMDL" << std::endl;
  return true;
}

}; // namespace trajan::io
