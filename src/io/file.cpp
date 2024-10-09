#include <algorithm>
#include <stdexcept>
#include <trajan/core/log.h>
#include <trajan/io/file.h>
#include <unordered_map>

namespace trajan::io {

using FileType = FileHandler::FileType;

void FileCompatibilityChecker::check_and_modify_handlers(
    std::vector<std::unique_ptr<FileHandler>> &handlers) {
  std::unordered_map<FileType, int> file_type_counts;

  for (const auto &handler : handlers) {
    file_type_counts[handler->get_file_type()]++;
  }
  if (file_type_counts.count(FileType::PDB) > 1) {
    trajan::log::warn("Multiple PDB files read. Using last.");
    auto it = std::find_if(handlers.begin(), handlers.end(),
                           [](const std::unique_ptr<FileHandler> &handler) {
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

std::unique_ptr<FileHandler> get_handler(std::string ext) {
  auto it = handler_map.find(ext);
  if (it != handler_map.end()) {
    std::unique_ptr<io::FileHandler> handler = it->second();
    trajan::log::debug(fmt::format("Recognised '{}' file extension", ext));
    return handler;
  } else {
    throw std::runtime_error(fmt::format("Unknown file extension: '{}'", ext));
  }
};

} // namespace trajan::io
