#pragma once
#include "trajan/io/file.h"

namespace trajan::io {

class PDBHandler : public FileHandler {
public:
  // ~PDBHandler() override = default;
  // PDBHandler();
  FileType get_file_type() const override { return FileType::PDB; }
  bool read_frame(core::Frame &frame) override;

private:
  constexpr std::string_view PDB_FMT =
      "%6s%5d %4s%c%3s %c%4d%c   %8lf%8lf%8lf%6lf%6lf          %2s%2s";
};

}; // namespace trajan::io
