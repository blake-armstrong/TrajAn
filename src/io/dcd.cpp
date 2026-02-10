// #include "trajan/core/unit_cell.h"
#include "trajan/core/util.h"
#include <occ/core/linear_algebra.h>
#include <occ/core/units.h>
#include <occ/crystal/unitcell.h>
#include <trajan/io/dcd.h>

namespace trajan::io {

using occ::Vec3;
using occ::crystal::UnitCell;
using trajan::core::Frame;

bool DCDHandler::_initialise() {
  if (m_mode == Mode::Read) {
    m_infile.open(this->file_path(), std::ios::binary);
    if (!m_infile.is_open()) {
      return false;
    }
    return this->parse_dcd_header();
  } else {
    m_outfile.open(this->file_path(), std::ios::binary);
    return m_outfile.is_open();
  }
}

void DCDHandler::_finalise() {
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

bool DCDHandler::read_next_frame(Frame &frame) {
  if (m_current_frame >= m_total_frames) {
    return false;
  }

  bool success = this->parse_dcd(frame);
  if (success) {
    m_current_frame++;
  }

  return success;
}

bool DCDHandler::write_next_frame(const Frame &frame) {
  if (m_current_frame == 0) {
    if (!write_dcd_header(frame)) {
      return false;
    }
  }

  if (!write_dcd_frame(frame)) {
    return false;
  }

  std::streampos current_pos = m_outfile.tellp();

  m_outfile.seekp(8, std::ios::beg);

  if (!write_binary(static_cast<int32_t>(m_current_frame + 1))) {
    return false;
  }

  m_outfile.seekp(current_pos);

  m_outfile.flush();

  m_current_frame++;

  return true;
}

template <typename T> bool DCDHandler::read_binary(T &value) {
  m_infile.read(reinterpret_cast<char *>(&value), sizeof(T));
  return m_infile.good();
}

template <typename T> bool DCDHandler::write_binary(const T &value) {
  m_outfile.write(reinterpret_cast<const char *>(&value), sizeof(T));
  return m_outfile.good();
}

bool DCDHandler::read_fortran_record(std::vector<char> &buffer) {
  int32_t record_length;
  if (!read_binary(record_length)) {
    return false;
  }

  buffer.resize(record_length);
  m_infile.read(buffer.data(), record_length);
  if (!m_infile.good()) {
    return false;
  }

  int32_t trailing_length;
  if (!read_binary(trailing_length)) {
    return false;
  }

  return (record_length == trailing_length);
}

bool DCDHandler::write_fortran_record(const std::vector<char> &buffer) {
  int32_t record_length = static_cast<int32_t>(buffer.size());

  if (!write_binary(record_length)) {
    return false;
  }

  m_outfile.write(buffer.data(), record_length);
  if (!m_outfile.good()) {
    return false;
  }

  if (!write_binary(record_length)) {
    return false;
  }

  return true;
}

bool DCDHandler::skip_fortran_record() {
  int32_t record_length;
  if (!read_binary(record_length)) {
    return false;
  }

  m_infile.seekg(record_length, std::ios::cur);
  if (!m_infile.good()) {
    return false;
  }

  int32_t trailing_length;
  return read_binary(trailing_length) && (record_length == trailing_length);
}

bool DCDHandler::_parse_dcd_header() {
  // First record: Main header (84 bytes)
  std::vector<char> header_buffer;
  if (!read_fortran_record(header_buffer)) {
    throw std::runtime_error("Failed to read DCD header record");
  }

  if (header_buffer.size() < 84) {
    throw std::runtime_error("DCD header record too small");
  }

  trajan::log::debug("Header record size: {} bytes", header_buffer.size());

  char signature[4];
  std::memcpy(signature, header_buffer.data(), 4);
  trajan::log::debug("File signature: '{}'", signature);

  if (std::strncmp(signature, "CORD", 4) != 0) {
    throw std::runtime_error("Invalid DCD signature");
  }

  int32_t *int_data = reinterpret_cast<int32_t *>(header_buffer.data());

  m_total_frames = static_cast<size_t>(int_data[1]); // NFILE
  int32_t first_timestep = int_data[2];              // ISTART
  int32_t save_frequency = int_data[3];              // NSAVC
  int32_t total_timesteps = int_data[4];             // NSTEP
  int32_t num_fixed_atoms = int_data[9];             // NAMNF
  int32_t charmm_version = int_data[20];             // CHARMM version

  trajan::log::debug("Trajectory parameters:");
  trajan::log::debug("  Total frames:        {}", m_total_frames);
  trajan::log::debug("  First timestep:      {}", first_timestep);
  trajan::log::debug("  Save frequency:      {} (every {} steps)",
                     save_frequency, save_frequency);
  trajan::log::debug("  Total timesteps:     {}", total_timesteps);
  trajan::log::debug("  Fixed atoms:         {}", num_fixed_atoms);

  if (charmm_version != 0) {
    m_is_charmm_format = true;
    m_has_extra_block = (int_data[11] != 0);
    m_has_4d_coords = (int_data[12] == 1);

    float timestep_float;
    std::memcpy(&timestep_float, &int_data[10], sizeof(float));
    m_timestep = static_cast<double>(timestep_float);

    trajan::log::debug("File format: CHARMM (version {})", charmm_version);
    trajan::log::debug("  Unit cell data:      {}",
                       m_has_extra_block ? "present" : "absent");
    trajan::log::debug("  4D coordinates:      {}",
                       m_has_4d_coords ? "yes" : "no");
    trajan::log::debug("  Timestep:            {} (float)", m_timestep);
  } else {
    // X-PLOR format
    m_is_charmm_format = false;
    m_has_extra_block = false;
    m_has_4d_coords = false;

    double timestep_double;
    std::memcpy(&timestep_double, &int_data[10], sizeof(double));
    m_timestep = timestep_double;

    trajan::log::debug("File format: X-PLOR");
    trajan::log::debug("  Timestep:            {} (double)", m_timestep);
  }

  if (num_fixed_atoms < 0) {
    trajan::log::warn("Negative fixed atoms count: {} - file may be corrupted",
                      num_fixed_atoms);
  }
  if (num_fixed_atoms > 100000) {
    trajan::log::warn(
        "Unusually high fixed atoms count: {} - may indicate file corruption",
        num_fixed_atoms);
  }
  if (m_total_frames == 0) {
    trajan::log::warn(
        "Header indicates zero frames - trajectory may be incomplete");
  }
  if (m_total_frames > 1000000) {
    trajan::log::warn("Very large frame count: {} - ensure this is expected",
                      m_total_frames);
  }
  if (save_frequency <= 0) {
    trajan::log::warn("Invalid save frequency: {}", save_frequency);
  }

  std::vector<char> title_buffer;
  if (!read_fortran_record(title_buffer)) {
    throw std::runtime_error("Failed to read DCD title record");
  }

  if (title_buffer.size() >= 4) {
    int32_t num_titles = *reinterpret_cast<int32_t *>(title_buffer.data());
    trajan::log::debug("Title record:");
    trajan::log::debug("  Number of titles:    {}", num_titles);
    trajan::log::debug("  Record size:         {} bytes", title_buffer.size());

    if (num_titles > 0 && num_titles < 100) {
      for (int i = 0;
           i < num_titles && (4 + (i + 1) * 80) <= title_buffer.size(); ++i) {
        std::string title(
            title_buffer.data() + 4 + i * 80,
            std::min(size_t(80), title_buffer.size() - (4 + i * 80)));
        title.erase(title.find_last_not_of(" \0") + 1);
        if (!title.empty()) {
          trajan::log::debug("  Title[{}]:           '{}'", i, title);
        }
      }
    }

    if (num_titles < 0 || num_titles > 1000) {
      trajan::log::warn("Suspicious title count: {}", num_titles);
    }
  }

  std::vector<char> atom_buffer;
  // Allocate coordinate buffers
  if (!read_fortran_record(atom_buffer)) {
    throw std::runtime_error("Failed to read DCD atom count record");
  }

  if (atom_buffer.size() < 4) {
    throw std::runtime_error(fmt::format(
        "DCD atom count record too small: {} bytes", atom_buffer.size()));
  }

  m_num_atoms =
      static_cast<size_t>(*reinterpret_cast<int32_t *>(atom_buffer.data()));
  trajan::log::debug("System size:");
  trajan::log::debug("  Number of atoms:     {}", m_num_atoms);

  if (m_num_atoms == 0) {
    throw std::runtime_error("Zero atoms in system");
  }
  if (m_num_atoms > 10000000) {
    trajan::log::warn("Very large system: {} atoms - ensure this is expected",
                      m_num_atoms);
  }

  m_x_coords.resize(m_num_atoms);
  m_y_coords.resize(m_num_atoms);
  m_z_coords.resize(m_num_atoms);

  trajan::log::debug("Allocated coordinate buffers for {} atoms", m_num_atoms);

  return true;
}

bool DCDHandler::parse_dcd_header() {
  try {
    trajan::log::debug("Attempting to parse DCD header:");
    bool success = this->_parse_dcd_header();
    if (success) {
      trajan::log::debug("Successfully parsed DCD header");
    }
    return success;
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format("Exception while parsing DCD header: {}", e.what()));
  }
}

bool DCDHandler::_parse_dcd(core::Frame &frame) {
  if (m_is_charmm_format && m_has_extra_block) {
    std::vector<char> unitcell_buffer;
    if (!read_fortran_record(unitcell_buffer)) {
      throw std::runtime_error("Failed to read unit cell record");
    }

    // Parse unit cell data (6 doubles: a, b, c, alpha, beta, gamma)
    if (unitcell_buffer.size() >= 48) { // 6 * 8 bytes
      double *cell_data = reinterpret_cast<double *>(unitcell_buffer.data());

      double a = cell_data[0];
      double b = cell_data[2];
      double c = cell_data[5];
      double alpha = cell_data[4];
      double beta = cell_data[3];
      double gamma = cell_data[1];
      trajan::log::trace("Unit cell (raw): a={:.4f} b={:.4f} c={:.4f}", a, b,
                         c);
      trajan::log::trace("Unit cell (raw): α={:.4f} β={:.4f} γ={:.4f}", alpha,
                         beta, gamma);
      // Check if angles are stored as cosines (CHARMM/NAMD format)
      if (alpha >= -1.0 && alpha <= 1.0 && beta >= -1.0 && beta <= 1.0 &&
          gamma >= -1.0 && gamma <= 1.0) {
        // Angles stored as cosines - convert to radians
        alpha = std::acos(alpha);
        beta = std::acos(beta);
        gamma = std::acos(gamma);
        trajan::log::trace(
            "Unit cell angles interpreted as cosines, converted to radians");
      } else {
        alpha = occ::units::radians(alpha);
        beta = occ::units::radians(beta);
        gamma = occ::units::radians(gamma);
        trajan::log::trace(
            "Unit cell angles interpreted as degrees, converted to radians");
      }
      using occ::units::degrees;
      trajan::log::trace("Unit cell: a={:.4f} b={:.4f} c={:.4f} Angstroms", a,
                         b, c);
      trajan::log::trace("Unit cell: α={:.4f}({:.4f}) β={:.4f}({:.4f}) "
                         "γ={:.4f}({:.4f}) rad(deg)",
                         alpha, degrees(alpha), beta, degrees(beta), gamma,
                         degrees(gamma));
      if (trajan::util::unitcell_is_reasonable(a, b, c, alpha, beta, gamma)) {
        trajan::log::trace("Unit cell is reasonable");
        UnitCell uc = occ::crystal::triclinic_cell(a, b, c, alpha, beta, gamma);
        frame.set_unit_cell(uc);
      } else {
        trajan::log::trace("Unit cell is NOT reasonable. Not using it.");
      }
    }
  }

  // Read X coordinates
  std::vector<char> coord_buffer;
  if (!read_fortran_record(coord_buffer)) {
    throw std::runtime_error("Failed to read X coordinates");
  }

  if (coord_buffer.size() != m_num_atoms * sizeof(float)) {
    throw std::runtime_error("X coordinate record size mismatch");
  }
  std::memcpy(m_x_coords.data(), coord_buffer.data(), coord_buffer.size());

  // Read Y coordinates
  if (!read_fortran_record(coord_buffer)) {
    throw std::runtime_error("Failed to read Y coordinates");
  }

  if (coord_buffer.size() != m_num_atoms * sizeof(float)) {
    throw std::runtime_error("Y coordinate record size mismatch");
  }
  std::memcpy(m_y_coords.data(), coord_buffer.data(), coord_buffer.size());

  // Read Z coordinates
  if (!read_fortran_record(coord_buffer)) {
    throw std::runtime_error("Failed to read Z coordinates");
  }

  if (coord_buffer.size() != m_num_atoms * sizeof(float)) {
    throw std::runtime_error("Z coordinate record size mismatch");
  }
  std::memcpy(m_z_coords.data(), coord_buffer.data(), coord_buffer.size());

  // Skip 4th dimension if present
  if (m_is_charmm_format && m_has_4d_coords) {
    if (!skip_fortran_record()) {
      throw std::runtime_error("Failed to skip 4D coordinates");
    }
  }

  // Update frame with coordinates
  for (size_t i = 0; i < m_num_atoms; ++i) {
    Vec3 pos = {m_x_coords[i], m_y_coords[i], m_z_coords[i]};
    frame.update_atom_position(i, pos);
  }

  return true;
}

bool DCDHandler::parse_dcd(Frame &frame) {
  try {
    trajan::log::trace("Parsing next DCD frame");
    auto success = this->_parse_dcd(frame);
    if (success) {
      trajan::log::trace("Successfully parsed DCD frame");
    }
    return success;
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format("Exception while parsing DCD frame: {}", e.what()));
    return false;
  }
}

bool DCDHandler::write_dcd_header(const Frame &frame) {
  m_num_atoms = frame.num_atoms();

  // Create 84-byte main header
  std::vector<char> header_buffer(84, 0);
  int32_t *int_data = reinterpret_cast<int32_t *>(header_buffer.data());

  // Set signature
  std::memcpy(&int_data[0], "CORD", 4);

  // Set header fields
  int_data[1] = 0; // NFILE (will be updated as frames are written)
  int_data[2] = 0; // ISTART
  int_data[3] = 1; // NSAVC
  int_data[4] = 0; // NSTEP (total steps)

  // Initialize remaining fields to 0
  for (int i = 5; i <= 19; ++i) {
    int_data[i] = 0;
  }

  // Set CHARMM version to indicate CHARMM format
  int_data[20] = 24; // CHARMM version

  // Timestep
  int_data[10] = static_cast<int32_t>(1.0f);

  // Set unit cell flag
  int_data[11] = frame.has_unit_cell() ? 1 : 0; // Has unit cell data

  if (!write_fortran_record(header_buffer)) {
    return false;
  }

  // Write title record
  std::vector<char> title_buffer(84, 0);
  int32_t num_titles = 1;
  std::memcpy(title_buffer.data(), &num_titles, sizeof(int32_t));
  std::string title = "Created by TrajAn";
  std::memcpy(title_buffer.data() + 4, title.c_str(),
              std::min(title.size(), size_t(79)));

  if (!write_fortran_record(title_buffer)) {
    return false;
  }

  // Write atom count record
  std::vector<char> atom_count_buffer(4);
  int32_t num_atoms = static_cast<int32_t>(m_num_atoms);
  std::memcpy(atom_count_buffer.data(), &num_atoms, sizeof(int32_t));

  return write_fortran_record(atom_count_buffer);
}

bool DCDHandler::write_unit_cell_data(const Frame &frame) {
  const auto &uc = frame.unit_cell();

  std::vector<char> unit_cell_buffer(48, 0); // 6 doubles = 48 bytes
  double *cell_data = reinterpret_cast<double *>(unit_cell_buffer.data());

  // Store unit cell parameters
  if (uc.has_value()) {
    cell_data[0] = uc->a();               // a
    cell_data[1] = std::cos(uc->gamma()); // cos(gamma) - angle between A and B
    cell_data[2] = uc->b();               // b
    cell_data[3] = std::cos(uc->beta());  // cos(beta) - angle between A and C
    cell_data[4] = std::cos(uc->alpha()); // cos(alpha) - angle between B and C
    cell_data[5] = uc->c();               // c
  } else {
    cell_data[0] = 0.0;
    cell_data[1] = 0.0;
    cell_data[2] = 0.0;
    cell_data[3] = 0.0;
    cell_data[4] = 0.0;
    cell_data[5] = 0.0;
  }

  return write_fortran_record(unit_cell_buffer);
}

bool DCDHandler::write_dcd_frame(const Frame &frame) {
  if (!write_unit_cell_data(frame)) {
    return false;
  }

  const auto &atoms = frame.atoms();
  m_x_coords.resize(m_num_atoms);
  m_y_coords.resize(m_num_atoms);
  m_z_coords.resize(m_num_atoms);

  for (size_t i = 0; i < m_num_atoms; ++i) {
    m_x_coords[i] = static_cast<float>(atoms[i].x);
    m_y_coords[i] = static_cast<float>(atoms[i].y);
    m_z_coords[i] = static_cast<float>(atoms[i].z);
  }

  std::vector<char> coord_buffer(m_num_atoms * sizeof(float));

  std::memcpy(coord_buffer.data(), m_x_coords.data(), coord_buffer.size());
  if (!write_fortran_record(coord_buffer)) {
    return false;
  }

  std::memcpy(coord_buffer.data(), m_y_coords.data(), coord_buffer.size());
  if (!write_fortran_record(coord_buffer)) {
    return false;
  }

  std::memcpy(coord_buffer.data(), m_z_coords.data(), coord_buffer.size());
  if (!write_fortran_record(coord_buffer)) {
    return false;
  }

  return true;
}

}; // namespace trajan::io
