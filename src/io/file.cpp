#include <algorithm>
#include <memory>
#include <stdexcept>
#include <trajan/core/log.h>
#include <trajan/io/file.h>
#include <unordered_map>

namespace trajan::io {

namespace fs = std::filesystem;
using FileType = FileHandler::FileType;
using uFileHandler = std::unique_ptr<FileHandler>;

void check_handlers(std::vector<uFileHandler> &handlers) {
  std::unordered_map<FileType, int> file_type_counts;

  for (const auto &handler : handlers) {
    file_type_counts[handler->get_file_type()]++;
  }
  if (file_type_counts.count(FileType::PDB) > 1) {
    trajan::log::warn("Multiple PDB files read. Using last.");
    auto it = std::find_if(handlers.begin(), handlers.end(),
                           [](const uFileHandler &handler) {
                             return handler->get_file_type() == FileType::PDB;
                           });
    if (it != handlers.end()) {
      handlers.erase(it);
      file_type_counts[FileType::PDB]--;
    }
  }
  if (file_type_counts.count(FileType::PDB) == 0 &&
      file_type_counts.count(FileType::DCD) > 0) {
    trajan::log::warn("DCD(s) loaded without a topology file.");
  }
  if (file_type_counts.count(FileType::PDB) == 0 &&
      file_type_counts.count(FileType::DCD) == 0) {
    throw std::runtime_error("No files!");
  }
};

uFileHandler get_handler(std::string ext) {
  auto it = handler_map.find(ext);
  if (it != handler_map.end()) {
    std::unique_ptr<io::FileHandler> handler = it->second();
    trajan::log::debug(fmt::format("Recognised '{}' file extension", ext));
    return handler;
  } else {
    throw std::runtime_error(fmt::format("Unknown file extension: '{}'", ext));
  }
};

uFileHandler read_input_file(const fs::path &file) {
  std::string ext = file.extension().string();
  trajan::log::debug("Attempting to read input from {}, file extension = {}",
                     file.generic_string(), ext);
  if (!fs::exists(file)) {
    throw std::runtime_error(
        fmt::format("Input file does not exist: '{}'", file.generic_string()));
  };

  uFileHandler handler = get_handler(ext);
  handler->set_file_path(file);

  return handler;
}

std::vector<uFileHandler> read_input_files(const std::vector<fs::path> &files) {
  std::vector<uFileHandler> handlers;
  for (const fs::path &file : files) {
    handlers.push_back(read_input_file(file));
  }
  check_handlers(handlers);
  return handlers;
}

} // namespace trajan::io
