#pragma once

#include <algorithm>
#include <optional>
#include <string>
#include <trajan/core/util.h>
#include <variant>
#include <vector>

namespace trajan::io {

struct IndexSelection {
  std::vector<int> indices;
};

struct AtomTypeSelection {
  std::vector<std::string> types;
};

struct MoleculeSelection {
  std::vector<int> molecule_ids;
};

using SelectionCriteria =
    std::variant<IndexSelection, AtomTypeSelection, MoleculeSelection>;

template <typename T> struct SelectionTraits {};

template <> struct SelectionTraits<IndexSelection> {
  using value_type = int;
  static constexpr bool allows_ranges = true;
  static constexpr char prefix = 'i';

  static std::optional<value_type> validate(const std::string &token) {
    try {
      int value = std::stoi(token);
      return value; // Allow any integer
    } catch (...) {
      return std::nullopt;
    }
  }
};

template <> struct SelectionTraits<MoleculeSelection> {
  using value_type = int;
  static constexpr bool allows_ranges = true;
  static constexpr char prefix = 'm';

  static std::optional<value_type> validate(const std::string &token) {
    try {
      int value = std::stoi(token);
      return value >= 0 ? std::optional<value_type>(value) : std::nullopt;
    } catch (...) {
      return std::nullopt;
    }
  }
};

template <> struct SelectionTraits<AtomTypeSelection> {
  using value_type = std::string;
  static constexpr bool allows_ranges = false;
  static constexpr char prefix = 't';

  static std::optional<value_type> validate(const std::string &token) {
    if (std::all_of(token.begin(), token.end(), [](char c) {
          return std::isalnum(c) || c == '_' || c == '*';
        })) {
      return token;
    }
    return std::nullopt;
  }
};

class SelectionParser {
public:
  static std::optional<SelectionCriteria> parse(const std::string &input);

private:
  template <typename SelectionType>
  static std::optional<SelectionCriteria>
  parse_selection(const std::string &input);

  template <typename SelectionType>
  static std::optional<
      std::pair<typename SelectionTraits<SelectionType>::value_type,
                typename SelectionTraits<SelectionType>::value_type>>
  parse_range(const std::string &input) {
    auto parts = trajan::util::split_string(input, '-');
    if (parts.size() != 2)
      return std::nullopt;

    auto start = SelectionTraits<SelectionType>::validate(parts[0]);
    auto end = SelectionTraits<SelectionType>::validate(parts[1]);

    if (!start || !end || *start > *end)
      return std::nullopt;
    return std::make_pair(*start, *end);
  }
};

inline auto selection_validator =
    [](std::optional<SelectionCriteria> &parsed_sel) {
      return [parsed_sel = &parsed_sel](const std::string &input) {
        auto result = SelectionParser::parse(input);
        if (!result) {
          return std::string("Invalid selection format");
        }
        *parsed_sel = result;
        return std::string();
      };
    };

} // namespace trajan::io
