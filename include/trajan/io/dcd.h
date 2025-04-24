#pragma once
#include <trajan/io/file.h>

namespace trajan::io {

class DCDHandler : public FileHandler {
public:
  inline FileType file_type() const override { return FileType::DCD; }
  bool initialise() override;
  void finalise() override;
  bool read_next_frame(core::Frame &frame) override;

private:
  bool parse_dcd_header();
  bool parse_dcd(core::Frame &frame);
  size_t m_current_frame;
  size_t m_total_frames;
};

} // namespace trajan::io
