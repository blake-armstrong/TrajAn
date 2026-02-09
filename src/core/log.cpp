#include <iostream>
#include <memory>
#include <occ/core/util.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <trajan/core/log.h>
#include <unistd.h>

namespace trajan::log {

namespace {
std::shared_ptr<spdlog::logger> m_currentlogger = spdlog::default_logger();

spdlog::level::level_enum verbosity_to_level(const std::string &verbosity) {
  std::string level_lower = occ::util::to_lower_copy(verbosity);
  if (level_lower == "debug")
    return spdlog::level::trace;
  if (level_lower == "verbose")
    return spdlog::level::debug;
  if (level_lower == "minimal")
    return spdlog::level::warn;
  if (level_lower == "silent")
    return spdlog::level::critical;
  return spdlog::level::info; // default for "normal" and unknown values
}

spdlog::level::level_enum verbosity_to_level(int verbosity) {
  switch (verbosity) {
  case 4:
    return spdlog::level::trace;
  case 3:
    return spdlog::level::debug;
  case 1:
    return spdlog::level::warn;
  case 0:
    return spdlog::level::critical;
  default:
    return spdlog::level::info;
  }
}

std::string level_to_pattern(spdlog::level::level_enum &n) {
  switch (n) {
  case spdlog::level::trace:
    return PATTERN_VERBOSE;
  case spdlog::level::debug:
    return PATTERN_VERBOSE;
  case spdlog::level::warn:
    return PATTERN_SIMPLE;
  case spdlog::level::critical:
    return PATTERN_SIMPLE;
  default:
    return PATTERN_SIMPLE;
  }
}
} // namespace

void set_log_level(spdlog::level::level_enum level) {
  m_currentlogger->set_level(level);
  std::string pattern = level_to_pattern(level);
  spdlog::set_pattern(pattern);
  spdlog::enable_backtrace(32);
}

void set_log_level(const std::string &verbosity) {
  auto level = verbosity_to_level(verbosity);
  m_currentlogger->set_level(level);
  std::string pattern = level_to_pattern(level);
  spdlog::set_pattern(pattern);
  spdlog::enable_backtrace(32);
}

void set_log_level(int verbosity) {
  auto level = verbosity_to_level(verbosity);
  m_currentlogger->set_level(level);
  std::string pattern = level_to_pattern(level);
  spdlog::set_pattern(pattern);
  spdlog::enable_backtrace(32);
}

void set_log_file(const std::string &filename) {
  try {
    auto file_logger = spdlog::basic_logger_mt("trajan_logger", filename,
                                               true); // true = truncate
    m_currentlogger = file_logger;
    file_logger->set_level(m_currentlogger->level());
    spdlog::set_default_logger(m_currentlogger);
  } catch (const spdlog::spdlog_ex &ex) {
    spdlog::warn(
        "Failed to create file logger: {}. Using existing logger instead.",
        ex.what());
  }
  spdlog::set_pattern(PATTERN_SIMPLE);
  spdlog::enable_backtrace(32);
}

void set_subcommand_log_pattern(const std::string &subcommand) {
  auto level = m_currentlogger->level();
  if (level <= spdlog::level::debug) {
    spdlog::set_pattern(fmt::format("[%H:%M:%S]-[{}]-[%^%l%$] %v", subcommand));
  }
}

Progress::Progress(int total, const std::string &label)
    : m_label(label), m_current(0), m_total(total), m_last_length(0),
      m_finished(false) {

  auto level = m_currentlogger->level();
  bool is_tty = isatty(fileno(stdout)) != 0;
  m_enabled = is_tty && level <= spdlog::level::info;

  if (m_enabled) {
    display();
  }
}

Progress::Progress(const std::string &label)
    : m_label(label), m_current(0), m_total(-1), m_last_length(0),
      m_finished(false) {

  auto level = m_currentlogger->level();
  bool is_tty = isatty(fileno(stdout)) != 0;
  m_enabled = is_tty && level <= spdlog::level::info;
}

Progress::~Progress() {
  if (!m_finished) {
    finish();
  }
}

void Progress::update(int current) {
  m_current = current;
  if (m_enabled && !m_finished) {
    display();
  }
}

void Progress::update(int current, const std::string &message) {
  m_current = current;
  if (m_enabled && !m_finished) {
    display_with_message(message);
  }
}

void Progress::increment() { update(m_current + 1); }

void Progress::finish() {
  if (m_enabled && !m_finished) {
    clear_line();
    m_finished = true;
  }
}

void Progress::finish(const std::string &final_message) {
  if (m_enabled && !m_finished) {
    clear_line();
    spdlog::info(final_message);
    m_finished = true;
  }
}

void Progress::clear_line() {
  if (m_last_length > 0) {
    std::cout << fmt::format("\r{}\r", std::string(m_last_length, ' '))
              << std::flush;
    m_last_length = 0;
  }
}

void Progress::display() {
  clear_line();

  std::string output;

  if (m_total > 0) {
    double percentage = (100.0 * m_current / m_total);
    const int bar_width = 40;
    int filled = static_cast<int>(bar_width * m_current / m_total);

    std::string bar;
    bar.reserve(bar_width);
    for (int i = 0; i < bar_width; ++i) {
      if (i < filled)
        bar += "=";
      else if (i == filled)
        bar += ">";
      else
        bar += " ";
    }

    output = fmt::format("\r{}: {:6d}/{} [{}] {:6.2f}%", m_label, m_current,
                         m_total, bar, percentage);
  } else {
    output = fmt::format("\r{}: {}", m_label, m_current);
  }

  m_last_length = output.length();
  std::cout << output << std::flush;
}

void Progress::display_with_message(const std::string &message) {
  clear_line();

  std::string output;

  if (m_total > 0) {
    double percentage = (100.0 * m_current / m_total);
    output = fmt::format("\r{}: {:6d}/{} [{:6.2f}%] {}", m_label, m_current,
                         m_total, percentage, message);
  } else {
    output = fmt::format("\r{}: {} {}", m_label, m_current, message);
  }
  m_last_length = output.length();
  std::cout << output << std::flush;
}
} // namespace trajan::log
