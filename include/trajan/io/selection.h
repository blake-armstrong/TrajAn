#pragma once

#include "spdlog/common.h"
#include <algorithm>
#include <optional>
#include <string>
#include <trajan/core/atom.h>
#include <trajan/core/log.h>
#include <trajan/core/molecule.h>
#include <trajan/core/neigh.h>
#include <trajan/core/util.h>
#include <variant>
#include <vector>

namespace trajan::io {

using Atom = trajan::core::Atom;
using Molecule = trajan::core::Molecule;

struct SelectionBase {
  virtual std::string name() const = 0;
  virtual ~SelectionBase() = default;
  std::string raw_input() const { return m_raw_input; }
  std::string m_raw_input{};
};

template <typename T> struct Selection : public SelectionBase {
  Selection(std::vector<T> d, const std::string &raw_input) {
    data = std::move(d);
    m_raw_input = raw_input;
  }
  std::vector<T> data;

  auto begin() { return data.begin(); }
  auto end() { return data.end(); }

  auto begin() const { return data.begin(); }
  auto end() const { return data.end(); }
};

struct AtomIndexSelection : public Selection<int> {
  using Selection<int>::Selection;
  std::string name() const override { return "AtomIndexSelection"; }
};

struct AtomTypeSelection : public Selection<std::string> {
  using Selection<std::string>::Selection;
  std::string name() const override { return "AtomTypeSelection"; }
};

struct MoleculeIndexSelection : public Selection<int> {
  using Selection<int>::Selection;
  std::string name() const override { return "MoleculeIndexSelection"; }
};

struct MoleculeTypeSelection : public Selection<std::string> {
  using Selection<std::string>::Selection;
  std::string name() const override { return "MoleculeTypeSelection"; }
};

using SelectionCriteria =
    std::variant<AtomIndexSelection, AtomTypeSelection, MoleculeIndexSelection,
                 MoleculeTypeSelection>;

template <typename T> struct SelectionTraits {};

template <> struct SelectionTraits<AtomIndexSelection> {
  using value_type = int;
  static constexpr bool allows_ranges = true;
  static constexpr char prefix = 'i';

  static std::optional<value_type> validate(const std::string &token) {
    try {
      int value = std::stoi(token);
      return value;
    } catch (...) {
      return std::nullopt;
    }
  }
};

template <> struct SelectionTraits<MoleculeIndexSelection> {
  using value_type = int;
  static constexpr bool allows_ranges = true;
  static constexpr char prefix = 'j';

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
  static constexpr char prefix = 'a';

  static std::optional<value_type> validate(const std::string &token) {
    if (std::all_of(token.begin(), token.end(), [](char c) {
          return std::isalnum(c) || c == '_' || c == '*';
        })) {
      return token;
    }
    return std::nullopt;
  }
};

template <> struct SelectionTraits<MoleculeTypeSelection> {
  using value_type = std::string;
  static constexpr bool allows_ranges = false;
  static constexpr char prefix = 'm';

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
  static std::optional<std::vector<SelectionCriteria>>
  parse(const std::string &input);

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

void process_selection(const SelectionCriteria &selection,
                       const std::vector<Atom> &atoms,
                       const std::vector<Molecule> &molecules,
                       std::vector<core::EntityVariant> &entities);

void process_selection(const std::vector<SelectionCriteria> &selections,
                       const std::vector<Atom> &atoms,
                       const std::vector<Molecule> &molecules,
                       std::vector<core::EntityVariant> &entities);

void update_entities_with_positions(std::vector<core::EntityVariant> &entities,
                                    const std::vector<Atom> &atoms,
                                    const std::vector<Molecule> &molecules);

void print_parsed_selection(const std::vector<SelectionCriteria> &result);

std::vector<SelectionCriteria> selection_validator(
    const std::string &input,
    std::optional<std::vector<char>> restrictions = std::nullopt);

} // namespace trajan::io
