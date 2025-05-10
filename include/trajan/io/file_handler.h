#pragma once

#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <trajan/core/frame.h>
#include <trajan/core/trajectory.h>
#include <vector>

namespace trajan::io {

namespace core = trajan::core;
namespace fs = std::filesystem;

class FileHandler {
public:
  enum class FileType {
    PDB,
    DCD,
  };
  virtual ~FileHandler() = default;

  inline virtual FileType file_type() const = 0;

  bool initialise();
  bool read_frame(core::Frame &frame);
  void finalise();

  // bool read_frame(core::Frame &frame);

  void set_file_path(fs::path file_path) {
    m_file_path = file_path;
    m_file_name = file_path.generic_string();
  }
  fs::path file_path() { return m_file_path; }
  std::string file_name() { return m_file_name; }

protected:
  std::ifstream m_file;
  virtual bool _initialise() = 0;
  virtual void _finalise() = 0;
  virtual bool read_next_frame(core::Frame &frame) = 0;
  bool validate_frame(core::Frame &frame);

private:
  fs::path m_file_path;
  std::string m_file_name;
};

using FileHandlerPtr = std::unique_ptr<FileHandler>;
using FileType = FileHandler::FileType;

FileHandlerPtr get_handler(std::string ext);
void check_handlers(std::vector<FileHandlerPtr> &handlers);
FileHandlerPtr read_input_file(const fs::path &file);
std::vector<FileHandlerPtr>
read_input_files(const std::vector<fs::path> &files);
}; // namespace trajan::io

namespace std {
template <> struct hash<trajan::io::FileType> {
  std::size_t operator()(const trajan::io::FileType &k) const {
    return static_cast<std::size_t>(k);
  }
};

}; // namespace std
