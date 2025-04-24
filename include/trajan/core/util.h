#pragma once
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <sstream>
#include <string>
#include <trajan/core/linear_algebra.h>
#include <vector>

namespace trajan::util {

namespace fs = std::filesystem;

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

static inline double square_distance(Vec3 &p1, Vec3 &p2) {
  double dx = p1.x() - p2.x(), dy = p1.y() - p2.y(), dz = p1.z() - p2.z();
  return dx * dx + dy * dy + dz * dz;
}

template <typename T>
static inline const std::vector<T>
combine_vectors(const std::vector<std::vector<T>> &vecs) {
  std::vector<T> combined_vec;
  for (const std::vector<T> &vec : vecs) {
    combined_vec.reserve(combined_vec.size() + vec.size());
    combined_vec.insert(combined_vec.end(), vec.begin(), vec.end());
  }
  return combined_vec;
}

struct Opts {
  size_t num_threads;
  std::vector<fs::path> infiles;
};

} // namespace trajan::util
