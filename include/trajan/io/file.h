#pragma once

#include <filesystem>
#include <memory>
#include <string>
#include <trajan/core/trajectory.h>
#include <unordered_map>
#include <vector>

namespace trajan::io {

namespace core = trajan::core;

class FileHandler {
public:
  enum class FileType {
    PDB,
    DCD,
  };
  std::filesystem::path file_path;
  virtual ~FileHandler() = default;
  virtual FileType get_file_type() const = 0;

  virtual core::Trajectory read_trajectory(const std::string &filename) = 0;

  virtual void write_trajectory(const std::string &filename,
                                const core::Trajectory &trajectory) = 0;
};

class FileCompatibilityChecker {
public:
  void check_and_modify_handlers(
      std::vector<std::unique_ptr<FileHandler>> &handlers);
};

class PDBHandler : public FileHandler {
public:
  // ~PDBHandler() override = default;
  // PDBHandler();
  FileType get_file_type() const override { return FileType::PDB; }
  core::Trajectory read_trajectory(const std::string &filename) override;
  void write_trajectory(const std::string &filename,
                        const core::Trajectory &trajectory) override;
};

class DCDHandler : public FileHandler {
public:
  // DCDHandler() = default;
  // ~DCDHandler() override = default;
  FileType get_file_type() const override { return FileType::DCD; }
  core::Trajectory read_trajectory(const std::string &filename) override;
  void write_trajectory(const std::string &filename,
                        const core::Trajectory &trajectory) override;
};

static const std::unordered_map<
    std::string, std::function<std::unique_ptr<io::FileHandler>()>>
    handler_map = {{".pdb", []() { return std::make_unique<PDBHandler>(); }},
                   {".dcd", []() { return std::make_unique<DCDHandler>(); }}};

std::unique_ptr<FileHandler> get_handler(std::string ext); //{
//   auto it = handlerMap.find(ext);
//   if (it != handlerMap.end()) {
//       trajan::log::debug("Detected " + ext + " input file");
//       std::unique_ptr<io::FileHandler> handler = it->second();
//       handler->file_name = filename;
//       return handler;
//   } else {
//       throw std::runtime_error(fmt::format("Unknown file type: '{}'",
//       filename));
//   }
// }

}; // namespace trajan::io

namespace std {
template <> struct hash<trajan::io::FileHandler::FileType> {
  std::size_t operator()(const trajan::io::FileHandler::FileType &k) const {
    return static_cast<std::size_t>(k);
  }
};

}; // namespace std
