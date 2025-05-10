#include <trajan/io/dcd.h>

namespace trajan::io {

namespace core = trajan::core;

bool DCDHandler::_initialise() {
  m_file.open(this->file_path(), std::ios::binary);
  if (!m_file.is_open()) {
    return false;
  }
  return this->parse_dcd_header();
}

void DCDHandler::_finalise() {
  if (m_file.is_open()) {
    m_file.close();
  }
}

bool DCDHandler::read_next_frame(core::Frame &frame) {
  if (m_current_frame >= m_total_frames) {
    return false;
  }

  bool success = this->parse_dcd(frame);

  // Increment frame counter
  if (success) {
    m_current_frame++;
  }

  return success;
}

bool DCDHandler::parse_dcd_header() {
  // DCD-specific header parsing logic
  // Set m_total_frames based on header
  return true;
}

bool DCDHandler::parse_dcd(core::Frame &frame) {
  // TODO:
  return true;
}

}; // namespace trajan::io
