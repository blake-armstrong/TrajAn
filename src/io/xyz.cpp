#include <fmt/core.h>
#include <occ/core/linear_algebra.h>
#include <occ/crystal/unitcell.h>
#include <trajan/core/log.h>
#include <trajan/io/xyz.h>

#include <cmath>
#include <sstream>

namespace trajan::io {

using occ::Vec3;
using occ::crystal::triclinic_cell;
using occ::crystal::UnitCell;
using trajan::core::EnhancedAtom;
using trajan::core::Frame;

namespace {

// Parse the Lattice key from an extended XYZ comment line.
// Expects 9 space-separated values (3x3 matrix, rows = lattice vectors)
// enclosed in double-quotes: Lattice="ax ay az bx by bz cx cy cz"
std::optional<UnitCell> parse_lattice(const std::string &comment) {
  size_t lat_pos = comment.find("Lattice=");
  if (lat_pos == std::string::npos) {
    return std::nullopt;
  }

  size_t val_start = lat_pos + 8; // length of "Lattice="
  if (val_start >= comment.size()) {
    return std::nullopt;
  }

  std::string value;
  if (comment[val_start] == '"') {
    size_t quote_end = comment.find('"', val_start + 1);
    if (quote_end == std::string::npos) {
      trajan::log::warn("XYZ: malformed Lattice value (unclosed quote)");
      return std::nullopt;
    }
    value = comment.substr(val_start + 1, quote_end - val_start - 1);
  } else {
    size_t end = comment.find_first_of(" \t", val_start);
    value = comment.substr(val_start, end == std::string::npos ? end
                                                               : end - val_start);
  }

  std::istringstream iss(value);
  double vals[9];
  for (int i = 0; i < 9; i++) {
    if (!(iss >> vals[i])) {
      trajan::log::warn("XYZ: failed to parse 9 values from Lattice field");
      return std::nullopt;
    }
  }

  // Rows of the matrix are the three lattice vectors a, b, c
  Vec3 a_vec(vals[0], vals[1], vals[2]);
  Vec3 b_vec(vals[3], vals[4], vals[5]);
  Vec3 c_vec(vals[6], vals[7], vals[8]);

  double a = a_vec.norm();
  double b = b_vec.norm();
  double c = c_vec.norm();

  if (a < 1e-10 || b < 1e-10 || c < 1e-10) {
    trajan::log::warn("XYZ: zero-length lattice vector");
    return std::nullopt;
  }

  // Clamp dot products to [-1, 1] to avoid NaN from acos
  double cos_alpha = std::clamp(b_vec.dot(c_vec) / (b * c), -1.0, 1.0);
  double cos_beta  = std::clamp(a_vec.dot(c_vec) / (a * c), -1.0, 1.0);
  double cos_gamma = std::clamp(a_vec.dot(b_vec) / (a * b), -1.0, 1.0);

  return triclinic_cell(a, b, c,
                        std::acos(cos_alpha),
                        std::acos(cos_beta),
                        std::acos(cos_gamma));
}

} // anonymous namespace

bool XYZHandler::_initialise() {
  if (m_mode == Mode::Read) {
    m_infile.open(file_path());
    return m_infile.is_open();
  } else {
    m_outfile.open(file_path());
    return m_outfile.is_open();
  }
}

void XYZHandler::_finalise() {
  if (m_mode == Mode::Read) {
    if (m_infile.is_open()) {
      m_infile.close();
    }
  } else {
    if (m_outfile.is_open()) {
      m_outfile.flush();
      m_outfile.close();
    }
  }
}

bool XYZHandler::read_next_frame(Frame &frame) {
  std::string line;

  // Skip blank lines between frames
  while (std::getline(m_infile, line)) {
    if (line.find_first_not_of(" \t\r\n") != std::string::npos) {
      break;
    }
  }

  if (!m_infile || line.empty()) {
    return false;
  }

  // First line: atom count
  size_t num_atoms;
  try {
    num_atoms = std::stoul(line);
  } catch (const std::exception &) {
    trajan::log::warn("XYZ: failed to parse atom count from: '{}'", line);
    return false;
  }

  if (num_atoms == 0) {
    trajan::log::warn("XYZ: frame declares zero atoms");
    return false;
  }

  // Second line: comment (may carry extended XYZ key=value pairs)
  std::string comment;
  if (!std::getline(m_infile, comment)) {
    trajan::log::warn("XYZ: unexpected EOF reading comment line");
    return false;
  }

  auto uc = parse_lattice(comment);
  if (uc.has_value()) {
    frame.set_unit_cell(uc.value());
    trajan::log::debug("XYZ: parsed unit cell from Lattice key");
  }

  // Decide whether to create new atoms or just update positions
  bool create_atoms = (frame.num_atoms() == 0);
  if (!create_atoms && frame.num_atoms() != num_atoms) {
    trajan::log::warn("XYZ: atom count mismatch (frame has {}, file has {})",
                      frame.num_atoms(), num_atoms);
    return false;
  }

  std::vector<core::Atom> atoms;
  if (create_atoms) {
    atoms.resize(num_atoms);
  }

  for (size_t i = 0; i < num_atoms; ++i) {
    if (!std::getline(m_infile, line)) {
      trajan::log::warn("XYZ: unexpected EOF reading atom {}", i);
      return false;
    }

    std::istringstream iss(line);
    std::string element_sym;
    double x, y, z;

    if (!(iss >> element_sym >> x >> y >> z)) {
      trajan::log::warn("XYZ: failed to parse atom line: '{}'", line);
      return false;
    }

    if (create_atoms) {
      EnhancedAtom atom(element_sym, Vec3(x, y, z));
      atom.index = static_cast<int>(i);
      atom.uindex = static_cast<int>(i);
      atom.molecule_index = 0;
      atom.umolecule_index = 0;
      atom.type = element_sym;
      atoms[i] = atom;
    } else {
      Vec3 pos(x, y, z);
      frame.update_atom_position(i, pos);
    }
  }

  if (create_atoms) {
    frame.set_atoms(atoms);
  }

  trajan::log::debug("XYZ: read frame with {} atoms", num_atoms);
  return true;
}

bool XYZHandler::write_next_frame(const Frame &frame) {
  const auto &atoms = frame.atoms();
  const size_t num_atoms = atoms.size();

  if (num_atoms == 0) {
    trajan::log::warn("XYZ: no atoms to write");
    return false;
  }

  // Line 1: atom count
  m_outfile << num_atoms << "\n";

  // Line 2: comment with optional extended XYZ Lattice key
  const auto &uc = frame.unit_cell();
  if (uc.has_value()) {
    const auto &unit_cell = uc.value();
    // Derive lattice vectors from unit cell
    Vec3 a_vec = unit_cell.to_cartesian(Vec3(1.0, 0.0, 0.0));
    Vec3 b_vec = unit_cell.to_cartesian(Vec3(0.0, 1.0, 0.0));
    Vec3 c_vec = unit_cell.to_cartesian(Vec3(0.0, 0.0, 1.0));
    m_outfile << fmt::format(
        "Lattice=\"{:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} "
        "{:.6f}\" Properties=species:S:1:pos:R:3",
        a_vec[0], a_vec[1], a_vec[2], b_vec[0], b_vec[1], b_vec[2], c_vec[0],
        c_vec[1], c_vec[2]);
  } else {
    m_outfile << "Properties=species:S:1:pos:R:3";
  }
  m_outfile << "\n";

  // Atom lines: element x y z
  for (const auto &atom : atoms) {
    std::string sym = atom.element.symbol();
    if (sym == "Xx" || sym.empty()) {
      sym = atom.type.empty() ? "X" : atom.type;
    }
    m_outfile << fmt::format("{} {:.6f} {:.6f} {:.6f}\n", sym, atom.x, atom.y,
                             atom.z);
  }

  return m_outfile.good();
}

} // namespace trajan::io
