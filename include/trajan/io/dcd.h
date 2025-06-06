#pragma once
#include <trajan/io/file_handler.h>

namespace trajan::io {

class DCDHandler : public FileHandler {
public:
  inline FileType file_type() const override { return FileType::DCD; }
  bool read_next_frame(core::Frame &frame) override;

private:
  bool _initialise() override;
  void _finalise() override;
  bool _parse_dcd_header();
  bool parse_dcd_header();
  bool _parse_dcd(core::Frame &frame);
  bool parse_dcd(core::Frame &frame);

  size_t m_current_frame = 0;
  size_t m_total_frames = 0;
  size_t m_num_atoms = 0;

  template <typename T> bool read_binary(T &value);
  bool read_fortran_record(std::vector<char> &buffer);
  bool skip_fortran_record();

  int m_dcd_version = 0;
  bool m_has_extra_block = false;
  bool m_has_4d_coords = false;
  bool m_is_charmm_format = false;

  std::vector<float> m_x_coords;
  std::vector<float> m_y_coords;
  std::vector<float> m_z_coords;

  std::ifstream m_file;
};

} // namespace trajan::io
