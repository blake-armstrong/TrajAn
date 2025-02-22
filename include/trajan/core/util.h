#pragma once

#include <algorithm>
#include <cctype>
#include <sstream>
#include <string>
#include <trajan/core/linear_algebra.h>
#include <vector>

namespace trajan::util {

static inline void capitalize(std::string &s) {
  s[0] = std::toupper(s[0]);
  std::transform(s.begin() + 1, s.end(), s.begin() + 1,
                 [](unsigned char c) { return std::tolower(c); });
}

static inline std::string capitalize_copy(std::string s) {
  capitalize(s);
  return s;
}

static inline void to_upper(std::string &s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return std::toupper(c); });
}

static inline std::string to_upper_copy(std::string s) {
  to_upper(s);
  return s;
}

static inline void to_lower(std::string &s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return std::tolower(c); });
}

static inline std::string to_lower_copy(std::string s) {
  to_lower(s);
  return s;
}

// trim from start (in place)
static inline void ltrim(std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                                  [](int ch) { return !std::isspace(ch); }));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(),
                       [](int ch) { return !std::isspace(ch); })
              .base(),
          s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
  ltrim(s);
  rtrim(s);
}

// trim from start (copying)
static inline std::string ltrim_copy(std::string s) {
  ltrim(s);
  return s;
}

// trim from end (copying)
static inline std::string rtrim_copy(std::string s) {
  rtrim(s);
  return s;
}

// trim from both ends (copying)
static inline std::string trim_copy(std::string s) {
  trim(s);
  return s;
}

static inline std::vector<std::string> split_string(const std::string &input,
                                                    char delimiter) {
  std::vector<std::string> tokens;
  std::stringstream ss(input);
  std::string token;

  while (std::getline(ss, token, delimiter)) {
    token.erase(0, token.find_first_not_of(" \t"));
    token.erase(token.find_last_not_of(" \t") + 1);
    if (!token.empty()) {
      tokens.push_back(token);
    }
  }
  return tokens;
}

template <typename T>
constexpr bool is_close(T a, T b,
                        const T rtol = Eigen::NumTraits<T>::dummy_precision(),
                        const T atol = Eigen::NumTraits<T>::epsilon()) {
  return abs(a - b) <= (atol + rtol * abs(b));
}

} // namespace trajan::util
