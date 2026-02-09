#include <string>
#include <trajan/core/util.h>
#include <trajan/io/selection.h>

namespace trajan::io {

template <typename SelectionType>
std::optional<SelectionCriteria>
SelectionParser::parse_selection(const std::string &input) {
  using Traits = SelectionTraits<SelectionType>;
  using ValueType = typename Traits::value_type;

  try {
    std::vector<ValueType> values;
    auto tokens = trajan::util::split_string(input, ',');

    for (const auto &token : tokens) {
      if (token.find('-') != std::string::npos) {
        if (!Traits::allows_ranges)
          return std::nullopt;

        auto range = parse_range<SelectionType>(token);
        if (!range)
          return std::nullopt;

        if constexpr (std::is_integral_v<ValueType>) {
          for (ValueType i = range->first; i <= range->second; ++i) {
            values.push_back(i);
          }
        }
      } else {
        auto value = Traits::validate(token);
        if (!value)
          return std::nullopt;
        values.push_back(*value);
      }
    }

    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());

    if constexpr (std::is_same_v<SelectionType, AtomIndexSelection>) {
      return SelectionCriteria{AtomIndexSelection{std::move(values)}};
    } else if constexpr (std::is_same_v<SelectionType,
                                        MoleculeIndexSelection>) {
      return SelectionCriteria{MoleculeIndexSelection{std::move(values)}};
    } else if constexpr (std::is_same_v<SelectionType, AtomTypeSelection>) {
      return SelectionCriteria{AtomTypeSelection{std::move(values)}};
    } else if constexpr (std::is_same_v<SelectionType, MoleculeTypeSelection>) {
      return SelectionCriteria{MoleculeTypeSelection{std::move(values)}};
    }
  } catch (...) {
    return std::nullopt;
  }

  return std::nullopt;
}

std::optional<std::vector<SelectionCriteria>>
SelectionParser::parse(const std::string &input) {
  if (input.empty())
    return std::nullopt;

  std::vector<SelectionCriteria> results;
  std::string current_token;
  char current_prefix = '\0';

  for (size_t i = 0; i < input.length(); i++) {
    char c = input[i];

    // Check if this is a prefix character
    if ((c == SelectionTraits<AtomIndexSelection>::prefix ||
         c == SelectionTraits<AtomTypeSelection>::prefix ||
         c == SelectionTraits<MoleculeIndexSelection>::prefix ||
         c == SelectionTraits<MoleculeTypeSelection>::prefix) &&
        (i == 0 || input[i - 1] == ',')) {

      // Process previous token if exists
      if (!current_token.empty() && current_prefix != '\0') {
        std::optional<SelectionCriteria> result;
        switch (current_prefix) {
        case SelectionTraits<AtomIndexSelection>::prefix:
          result = parse_selection<AtomIndexSelection>(current_token);
          break;
        case SelectionTraits<AtomTypeSelection>::prefix:
          result = parse_selection<AtomTypeSelection>(current_token);
          break;
        case SelectionTraits<MoleculeIndexSelection>::prefix:
          result = parse_selection<MoleculeIndexSelection>(current_token);
          break;
        case SelectionTraits<MoleculeTypeSelection>::prefix:
          result = parse_selection<MoleculeTypeSelection>(current_token);
          break;
        }
        if (!result)
          return std::nullopt;
        results.push_back(*result);
        current_token.clear();
      }
      current_prefix = c;
    } else {
      current_token += c;
    }
  }

  // Process final token
  if (!current_token.empty() && current_prefix != '\0') {
    std::optional<SelectionCriteria> result;
    switch (current_prefix) {
    case SelectionTraits<AtomIndexSelection>::prefix:
      result = parse_selection<AtomIndexSelection>(current_token);
      break;
    case SelectionTraits<AtomTypeSelection>::prefix:
      result = parse_selection<AtomTypeSelection>(current_token);
      break;
    case SelectionTraits<MoleculeIndexSelection>::prefix:
      result = parse_selection<MoleculeIndexSelection>(current_token);
      break;
    }
    if (!result)
      return std::nullopt;
    results.push_back(*result);
  }

  return results.empty() ? std::nullopt : std::make_optional(results);
}

void print_parsed_selection(const std::vector<SelectionCriteria> &result) {
  std::vector<int> atom_indices;
  std::vector<std::string> atom_types;
  std::vector<int> molecule_indices;
  std::vector<std::string> molecule_types;

  for (const auto &sel : result) {
    std::visit(
        [&](const auto &s) {
          using SelType = std::decay_t<decltype(s)>;

          if constexpr (std::is_same_v<SelType, AtomIndexSelection>) {
            atom_indices.insert(atom_indices.end(), s.data.begin(),
                                s.data.end());
          } else if constexpr (std::is_same_v<SelType, AtomTypeSelection>) {
            atom_types.insert(atom_types.end(), s.data.begin(), s.data.end());
          } else if constexpr (std::is_same_v<SelType,
                                              MoleculeIndexSelection>) {
            molecule_indices.insert(molecule_indices.end(), s.data.begin(),
                                    s.data.end());
          } else if constexpr (std::is_same_v<SelType, MoleculeTypeSelection>) {
            molecule_types.insert(molecule_types.end(), s.data.begin(),
                                  s.data.end());
          }
        },
        sel);
  }

  trajan::log::debug("Selection Summary:");
  trajan::log::debug("  atom indices: {}",
                     trajan::util::format_vector(atom_indices, 60));
  trajan::log::debug("  atom types: {}",
                     trajan::util::format_vector(atom_types, 60));
  trajan::log::debug("  molecule indices: {}",
                     trajan::util::format_vector(molecule_indices, 60));
  trajan::log::debug("  molecule types: {}",
                     trajan::util::format_vector(molecule_types, 60));
};

std::vector<SelectionCriteria>
selection_validator(const std::string &input,
                    std::optional<std::vector<char>> restrictions) {

  trajan::log::debug("Raw selection string: '{}'", input);
  auto result = SelectionParser::parse(input);
  if (!result) {
    throw std::invalid_argument("Invalid selection format");
  }
  if (restrictions.has_value()) {
    for (const auto &sel : result.value()) {
      std::visit(
          [restrictions = &restrictions.value()](const auto &s) {
            using SelType = std::decay_t<decltype(s)>;
            const char p = SelectionTraits<SelType>::prefix;
            if (std::find(restrictions->begin(), restrictions->end(), p) ==
                restrictions->end()) {
              throw std::invalid_argument(
                  fmt::format("Selection prefix '{}' is not allowed.", p));
            }
          },
          sel);
    };
  }
  print_parsed_selection(result.value());
  return result.value();
};

} // namespace trajan::io
