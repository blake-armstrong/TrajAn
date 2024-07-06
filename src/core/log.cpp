#include <tacpp/core/log.h>
#include <tacpp/core/util.h>

namespace tacpp::log {

void setup_logging(const std::string &verbosity) {
  auto level = tacpp::log::level::info;
  std::string level_lower = tacpp::util::to_lower_copy(verbosity);
  if (level_lower == "debug")
    level = tacpp::log::level::trace;
  else if (level_lower == "normal")
    level = tacpp::log::level::info;
  else if (level_lower == "verbose")
    level = tacpp::log::level::debug;
  else if (level_lower == "minimal")
    level = tacpp::log::level::warn;
  else if (level_lower == "silent")
    level = tacpp::log::level::critical;
  tacpp::log::set_level(level);
  spdlog::set_level(level);
  // store the last 32 debug messages in a buffer
  spdlog::enable_backtrace(32);
  spdlog::set_pattern("%v");
}

void setup_logging(int verbosity) {
  auto level = tacpp::log::level::info;
  switch (verbosity) {
  case 4:
    level = tacpp::log::level::trace;
    break;
  case 3:
    level = tacpp::log::level::debug;
    break;
  case 1:
    level = tacpp::log::level::warn;
    break;
  case 0:
    level = tacpp::log::level::critical;
    break;
  default:
    level = tacpp::log::level::info;
    break;
  }
  tacpp::log::set_level(level);
  spdlog::set_level(level);
  // store the last 32 debug messages in a buffer
  spdlog::enable_backtrace(32);
  spdlog::set_pattern("%v");
}
} // namespace tacpp::log
