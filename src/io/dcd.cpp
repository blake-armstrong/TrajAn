#include <trajan/io/dcd.h>

namespace trajan::io {

namespace core = trajan::core;

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
      // Update the number of frames in the header
      // m_outfile.seekp(4, std::ios::beg);
      m_outfile.flush();
      m_outfile.close();
    }
  }
}

bool DCDHandler::read_next_frame(core::Frame &frame) {
  if (m_current_frame >= m_total_frames) {
    return false;
  }

  bool success = this->parse_dcd(frame);
  if (success) {
    m_current_frame++;
  }

  return success;
}

bool DCDHandler::write_next_frame(const core::Frame &frame) {
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
  int32_t record_length = buffer.size();
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
  std::vector<char> header_buffer;
  if (!read_fortran_record(header_buffer)) {
    throw std::runtime_error("Failed to read DCD header record");
    return false;
  }

  if (header_buffer.size() < 84) {
    throw std::runtime_error("DCD header record too small");
    return false;
  }

  char signature[4];
  std::memcpy(signature, header_buffer.data(), 4);
  if (std::strncmp(signature, "CORD", 4) != 0) {
    throw std::runtime_error("Invalid DCD signature");
    return false;
  }

  int32_t *int_data = reinterpret_cast<int32_t *>(header_buffer.data());

  m_total_frames = static_cast<size_t>(int_data[1]); // NFILE
  int32_t first_frame = int_data[2];                 // NPRIV
  int32_t frame_skip = int_data[3];                  // NSAVC
  int32_t total_steps = int_data[4];                 // NSTEP

  m_is_charmm_format = (int_data[11] != 0);
  m_has_extra_block = (int_data[11] & 0x01) != 0;
  m_has_4d_coords = (int_data[11] & 0x02) != 0;

  trajan::log::debug(fmt::format("DCD Header Info:\n"
                                 "  Total frames: {}\n"
                                 "  CHARMM format: {}\n",
                                 m_total_frames,
                                 m_is_charmm_format ? "Yes" : "No"));

  if (!read_fortran_record(header_buffer)) {
    throw std::runtime_error("Failed to read DCD title record");
    return false;
  }

  if (header_buffer.size() >= 4) {
    int32_t num_titles = *reinterpret_cast<int32_t *>(header_buffer.data());
    trajan::log::debug("Number of titles: {}", num_titles);
  }

  // Read third header record (number of atoms)
  if (!read_fortran_record(header_buffer)) {
    throw std::runtime_error("Failed to read DCD atom count record");
    return false;
  }

  if (header_buffer.size() < 4) {
    throw std::runtime_error("DCD atom count record too small");
    return false;
  }

  m_num_atoms =
      static_cast<size_t>(*reinterpret_cast<int32_t *>(header_buffer.data()));
  trajan::log::debug("Number of atoms: {}", m_num_atoms);

  m_x_coords.resize(m_num_atoms);
  m_y_coords.resize(m_num_atoms);
  m_z_coords.resize(m_num_atoms);

  return true;
}

bool DCDHandler::parse_dcd_header() {
  try {
    return this->_parse_dcd_header();
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format("Exception while parsing DCD header: {}", e.what()));
  }
}

bool DCDHandler::_parse_dcd(core::Frame &frame) {
  if (m_is_charmm_format && m_has_extra_block) {
    if (!skip_fortran_record()) {
      throw std::runtime_error("Failed to skip unit cell record");
      return false;
    }
  }

  std::vector<char> coord_buffer;
  if (!read_fortran_record(coord_buffer)) {
    throw std::runtime_error("Failed to read X coordinates");
    return false;
  }

  if (coord_buffer.size() != m_num_atoms * sizeof(float)) {
    throw std::runtime_error("X coordinate record size mismatch");
    return false;
  }

  std::memcpy(m_x_coords.data(), coord_buffer.data(), coord_buffer.size());

  if (!read_fortran_record(coord_buffer)) {
    throw std::runtime_error("Failed to read Y coordinates");
    return false;
  }

  if (coord_buffer.size() != m_num_atoms * sizeof(float)) {
    throw std::runtime_error("Y coordinate record size mismatch");
    return false;
  }

  std::memcpy(m_y_coords.data(), coord_buffer.data(), coord_buffer.size());

  if (!read_fortran_record(coord_buffer)) {
    throw std::runtime_error("Failed to read Z coordinates");
    return false;
  }

  if (coord_buffer.size() != m_num_atoms * sizeof(float)) {
    throw std::runtime_error("Z coordinate record size mismatch");
    return false;
  }

  std::memcpy(m_z_coords.data(), coord_buffer.data(), coord_buffer.size());

  if (m_has_4d_coords) {
    if (!skip_fortran_record()) {
      throw std::runtime_error("Failed to skip 4D coordinates");
      return false;
    }
  }

  for (size_t i = 0; i < m_num_atoms; ++i) {
    Vec3 pos = {m_x_coords[i], m_y_coords[i], m_z_coords[i]};
    frame.update_atom_position(i, pos);
  }

  return true;
}

bool DCDHandler::parse_dcd(core::Frame &frame) {
  try {
    return this->_parse_dcd(frame);
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format("Exception while parsing DCD frame: {}", e.what()));
    return false;
  }
}

bool DCDHandler::write_dcd_header(const core::Frame &frame) {
  m_num_atoms = frame.num_atoms();

  std::vector<char> header_buffer(84, 0);
  int32_t *int_data = reinterpret_cast<int32_t *>(header_buffer.data());

  // Set the signature as the first integer (4 bytes)
  std::memcpy(&int_data[0], "CORD", 4);

  int_data[1] = 0; // NFILE
  int_data[2] = 0; // NPRIV
  int_data[3] = 1; // NSAVC
  int_data[4] = 0; // NSTEP
  // Set remaining fields to 0
  for (int i = 5; i <= 10; ++i) {
    int_data[i] = 0;
  }
  int_data[11] = 0; // Set flags to 0 for now
  // Set remaining fields to 0
  for (int i = 12; i < 21; ++i) {
    int_data[i] = 0;
  }

  if (!write_fortran_record(header_buffer)) {
    return false;
  }

  std::vector<char> title_buffer(84, 0);
  int32_t num_titles = 1;
  std::memcpy(title_buffer.data(), &num_titles, sizeof(int32_t));
  std::string title = "Created by TrajAn";
  std::memcpy(title_buffer.data() + 4, title.c_str(), title.size());

  if (!write_fortran_record(title_buffer)) {
    return false;
  }

  std::vector<char> atom_count_buffer(4, 0);
  int32_t num_atoms = m_num_atoms;
  std::memcpy(atom_count_buffer.data(), &num_atoms, sizeof(int32_t));

  return write_fortran_record(atom_count_buffer);
}

bool DCDHandler::write_dcd_frame(const core::Frame &frame) {
  const auto &atoms = frame.atoms();
  m_x_coords.resize(m_num_atoms);
  m_y_coords.resize(m_num_atoms);
  m_z_coords.resize(m_num_atoms);

  for (size_t i = 0; i < m_num_atoms; ++i) {
    m_x_coords[i] = atoms[i].x;
    m_y_coords[i] = atoms[i].y;
    m_z_coords[i] = atoms[i].z;
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
