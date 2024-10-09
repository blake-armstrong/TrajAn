#include <trajan/core/log.h>
#include <trajan/core/util.h>

namespace trajan::log {

void setup_logging(const std::string &verbosity) {
  auto level = trajan::log::level::info;
  std::string level_lower = trajan::util::to_lower_copy(verbosity);
  if (level_lower == "debug")
    level = trajan::log::level::trace;
  else if (level_lower == "normal")
    level = trajan::log::level::info;
  else if (level_lower == "verbose")
    level = trajan::log::level::debug;
  else if (level_lower == "minimal")
    level = trajan::log::level::warn;
  else if (level_lower == "silent")
    level = trajan::log::level::critical;
  trajan::log::set_level(level);
  spdlog::set_level(level);
  // store the last 32 debug messages in a buffer
  spdlog::enable_backtrace(32);
  spdlog::set_pattern("%v");
}

void setup_logging(int verbosity) {
  auto level = trajan::log::level::info;
  switch (verbosity) {
  case 4:
    level = trajan::log::level::trace;
    break;
  case 3:
    level = trajan::log::level::debug;
    break;
  case 1:
    level = trajan::log::level::warn;
    break;
  case 0:
    level = trajan::log::level::critical;
    break;
  default:
    level = trajan::log::level::info;
    break;
  }
  trajan::log::set_level(level);
  spdlog::set_level(level);
  // store the last 32 debug messages in a buffer
  spdlog::enable_backtrace(32);
  spdlog::set_pattern("%v");
}
} // namespace trajan::log
