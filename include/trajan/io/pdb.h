#pragma once
#include <string_view>
#include <trajan/io/file_handler.h>

namespace trajan::io {

constexpr std::string_view PDB_CRYST_FMT_READ =
    "%6c%9lf%9lf%9lf%7lf%7lf%7lf%11c%4c";
constexpr std::string_view PDB_LINE_FMT_READ =
    // AAAAAABBBBBXCCCCDEEEXFGGGGHXXXIIIIIIIIJJJJJJJJKKKKKKKKLLLLLLMMMMMMXXXXXXXXXXNNOO
    "%6c%5d%1c%4c%1c%3c%1c%1c%4d%1c%3c%8lf%8lf%8lf%6f%6f%10c%2c%2c";
constexpr std::string_view PDB_CONECT_FMT_READ = "%6c%5d%5d%5d%5d%5d";

constexpr std::string_view PDB_LINE_FMT_WRITE =
    "{:<6}{:>5}{:1}{:>4}{:1}{:>3}{:1}{:1}{:>4}{:1}{:>3}{:8.3f}{:8.3f}{:8.3f}{:"
    "6.2f}{:6.2f}{:10}{:>2}{:>2}";

class PDBHandler : public FileHandler {
public:
  inline FileType file_type() const override { return FileType::PDB; }
  bool parse_pdb(core::Frame &frame);
  bool read_next_frame(core::Frame &frame) override;
  bool write_next_frame(const core::Frame &frame) override;
  void set_original_ids(bool v) override { m_original_ids = v; }

private:
  bool _initialise() override;
  void _finalise() override;
  // Cached unit cell: CRYST1 typically appears only in the first frame of a
  // multi-frame PDB, so we carry it forward to subsequent frames.
  std::optional<occ::crystal::UnitCell> m_cached_uc;
  // A line read ahead that belongs to the next frame (e.g. an ATOM record
  // encountered after ENDMDL, signalling a new frame with no MODEL separator).
  std::string m_pending_line{};
  bool m_has_pending{false};
  bool m_original_ids{false};
  int m_write_frame_count{0};
  // CONECT lines cached from the first frame that has topology. Written once
  // at the end of the file after all ENDMDL records.
  std::vector<std::string> m_conect_lines{};
};

} // namespace trajan::io
