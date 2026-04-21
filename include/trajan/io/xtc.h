#pragma once

#include <cstdint>
#include <fstream>
#include <vector>
#include <trajan/io/file_handler.h>

namespace trajan::io {

class XTCHandler : public FileHandler {
public:
  inline FileType file_type() const override { return FileType::XTC; }

protected:
  bool _initialise() override;
  void _finalise() override;
  bool read_next_frame(core::Frame &frame) override;
  bool write_next_frame(const core::Frame &frame) override;

private:
  std::ifstream  m_infile;
  std::ofstream  m_outfile;
  std::streampos m_file_size{0};
  uint32_t       m_natoms{0};
  uint32_t       m_step{0};
  float          m_precision{1000.0f};   // quantisation precision (1/nm)
  std::vector<float>   m_coords;         // interleaved xyz, nm
  std::vector<uint8_t> m_buf;            // compressed data buffer
  std::vector<int32_t> m_int_coords;     // quantised coords for compression
};

} // namespace trajan::io
