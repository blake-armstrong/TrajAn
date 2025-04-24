#pragma once
#include <string_view>
#include <trajan/io/file.h>

namespace trajan::io {

constexpr std::string_view PDB_CRYST_FMT = "%6s%9lf%9lf%9lf%7lf%7lf%7lf%11s%4s";
constexpr std::string_view PDB_LINE_FMT =
    "%6s%5d %4s%c%3s %c%4d%c   %8lf%8lf%8lf%6lf%6lf          %2s%2s";

class PDBHandler : public FileHandler {
public:
  inline FileType file_type() const override { return FileType::PDB; }
  bool initialise() override;
  void finalise() override;
  bool parse_pdb(core::Frame &frame);
  bool read_next_frame(core::Frame &frame) override;

private:
  bool m_has_read{false};
};

} // namespace trajan::io
