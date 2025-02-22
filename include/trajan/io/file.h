#pragma once

#include <filesystem>
#include <memory>
#include <string>
#include <trajan/core/frame.h>
#include <trajan/core/trajectory.h>
#include <unordered_map>
#include <vector>

namespace trajan::io {

namespace core = trajan::core;
namespace fs = std::filesystem;
constexpr std::string_view PDB_CRYST_FMT = "%6s%9lf%9lf%9lf%7lf%7lf%7lf%11s%4s";
constexpr std::string_view PDB_LINE_FMT =
    "%6s%5d %4s%c%3s %c%4d%c   %8lf%8lf%8lf%6lf%6lf          %2s%2s";

class FileHandler {
public:
  enum class FileType {
    PDB,
    DCD,
  };
  virtual ~FileHandler() = default;
  virtual FileType get_file_type() const = 0;

  virtual bool read_frame(core::Frame &frame) = 0;
  void set_file_path(fs::path file_path) {
    m_file_path = file_path;
    m_file_name = file_path.generic_string();
  }
  fs::path get_file_path() { return m_file_path; }
  std::string get_file_name() { return m_file_name; }

private:
  fs::path m_file_path;
  std::string m_file_name;
};

using uFileHandler = std::unique_ptr<FileHandler>;

class PDBHandler : public FileHandler {
public:
  FileType get_file_type() const override { return FileType::PDB; }
  bool read_frame(core::Frame &frame) override;
};

class DCDHandler : public FileHandler {
public:
  FileType get_file_type() const override { return FileType::DCD; }
  bool read_frame(core::Frame &frame) override;
};

static const std::unordered_map<std::string, std::function<uFileHandler()>>
    handler_map = {{".pdb", []() { return std::make_unique<PDBHandler>(); }},
                   {".dcd", []() { return std::make_unique<DCDHandler>(); }}};

uFileHandler get_handler(std::string ext);
void check_handlers(std::vector<uFileHandler> &handlers);
uFileHandler read_input_file(const fs::path &file);
std::vector<uFileHandler> read_input_files(const std::vector<fs::path> &files);
}; // namespace trajan::io

namespace std {
template <> struct hash<trajan::io::FileHandler::FileType> {
  std::size_t operator()(const trajan::io::FileHandler::FileType &k) const {
    return static_cast<std::size_t>(k);
  }
};

}; // namespace std
