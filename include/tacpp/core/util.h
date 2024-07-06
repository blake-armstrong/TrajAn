#pragma once

namespace tacpp::util {

static inline void to_lower(std::string &s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return std::tolower(c); });
}

static inline std::string to_lower_copy(std::string s) {
  to_lower(s);
  return s;
}

} // namespace tacpp::util
