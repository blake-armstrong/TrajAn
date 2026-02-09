#pragma once
#include <spdlog/spdlog.h>
#include <string>

namespace trajan::log {
using spdlog::critical;
using spdlog::debug;
using spdlog::error;
using spdlog::info;
using spdlog::set_level;
using spdlog::trace;
using spdlog::warn;

constexpr std::string PATTERN_VERBOSE = "[%H:%M:%S]-[%^%l%$] %v";
constexpr std::string PATTERN_SIMPLE = "%v";
constexpr size_t LINE_WIDTH = 100;

namespace level {
using spdlog::level::critical;
using spdlog::level::debug;
using spdlog::level::err;
using spdlog::level::info;
using spdlog::level::trace;
using spdlog::level::warn;
} // namespace level

void set_log_level(const std::string &verbosity);
void set_log_level(spdlog::level::level_enum level);
void set_log_level(int verbosity);
void set_subcommand_log_pattern(const std::string &subcommand);

void set_log_file(const std::string &filename);

inline void flush() { spdlog::default_logger()->flush(); }

inline void flush_on(spdlog::level::level_enum level) {
  spdlog::flush_on(level);
}

inline void flush_every(std::chrono::seconds interval) {
  spdlog::flush_every(interval);
}

class Progress {
public:
  Progress(int total, const std::string &label = "Progress");
  explicit Progress(const std::string &label = "Counter");
  ~Progress();
  Progress(const Progress &) = delete;
  Progress &operator=(const Progress &) = delete;
  Progress(Progress &&) noexcept = default;
  Progress &operator=(Progress &&) noexcept = default;

  void update(int current);
  void update(int current, const std::string &message);
  void increment();
  void finish();
  void finish(const std::string &final_message);

private:
  void clear_line();
  void display();
  void display_with_message(const std::string &message);

  std::string m_label;
  int m_current;
  int m_total; // -1 means no total (counter mode)
  int m_last_length;
  bool m_finished;
  bool m_enabled;
};

} // namespace trajan::log
