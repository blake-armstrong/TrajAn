#include "trajan/core/neigh.h"
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
      return SelectionCriteria{AtomIndexSelection{std::move(values), input}};
    } else if constexpr (std::is_same_v<SelectionType,
                                        MoleculeIndexSelection>) {
      return SelectionCriteria{
          MoleculeIndexSelection{std::move(values), input}};
    } else if constexpr (std::is_same_v<SelectionType, AtomTypeSelection>) {
      return SelectionCriteria{AtomTypeSelection{std::move(values), input}};
    } else if constexpr (std::is_same_v<SelectionType, MoleculeTypeSelection>) {
      return SelectionCriteria{MoleculeTypeSelection{std::move(values), input}};
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
    case SelectionTraits<MoleculeTypeSelection>::prefix:
      result = parse_selection<MoleculeTypeSelection>(current_token);
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

void process_selection(const SelectionCriteria &selection,
                       const std::vector<Atom> &atoms,
                       const std::vector<Molecule> &molecules,
                       std::vector<core::EntityVariant> &entities) {
  std::visit(
      [&](const auto &sel) {
        using SelType = std::decay_t<decltype(sel)>;
        trajan::log::debug("Processing selection '{}'", sel.raw_input());
        trajan::log::debug("Selection identified as '{}'", sel.name());

        if constexpr (std::is_same_v<SelType, io::AtomIndexSelection>) {
          trajan::log::debug(
              "The following atoms have been idenitified from {}:", sel.name());
          trajan::log::debug(" NB: atoms are only printed if log level is 4.");
          for (const Atom &atom : atoms) {
            for (const int &ai : sel) {
              if (atom.index == ai) {
                entities.push_back(atom);
                trajan::log::trace(" {}", atom.repr());
              }
            }
          }
        } else if constexpr (std::is_same_v<SelType, io::AtomTypeSelection>) {
          trajan::log::debug(
              "The following atoms have been idenitified from {}:", sel.name());
          trajan::log::debug(" NB: atoms are only printed if log level is 4.");
          for (const Atom &atom : atoms) {
            for (const std::string &at : sel) {
              if (atom.type == at) {
                entities.push_back(atom);
                trajan::log::trace(" {}", atom.repr());
              }
            }
          }
        } else if constexpr (std::is_same_v<SelType,
                                            io::MoleculeIndexSelection>) {
          trajan::log::debug(
              "The following molecules have been idenitified from {}:",
              sel.name());
          trajan::log::debug(" NB: atoms are only printed if log level is 4.");
          for (const Molecule &molecule : molecules) {
            for (const int &mi : sel) {
              if (molecule.index == mi) {
                entities.push_back(molecule);
                trajan::log::trace("  {}", molecule.repr());
              }
            }
          }
        } else if constexpr (std::is_same_v<SelType,
                                            io::MoleculeTypeSelection>) {
          trajan::log::debug(
              "The following molecules have been idenitified from {}:",
              sel.name());
          trajan::log::debug(" NB: atoms are only printed if log level is 4.");
          for (const Molecule &molecule : molecules) {
            for (const std::string &mt : sel) {
              if (molecule.type == mt) {
                entities.push_back(molecule);
                trajan::log::trace("  {}", molecule.repr());
              }
            }
          }
        }
        size_t num_entities = entities.size();
        if (num_entities == 0) {
          throw std::runtime_error(fmt::format(
              "No entities found in selection '{}'", sel.raw_input()));
        };
        trajan::log::debug("Identified {} entities in selection '{}'",
                           num_entities, sel.raw_input());
      },
      selection);
};

void process_selection(const std::vector<SelectionCriteria> &selections,
                       const std::vector<Atom> &atoms,
                       const std::vector<Molecule> &molecules,
                       std::vector<core::EntityVariant> &entities) {

  for (const auto &selection : selections) {
    process_selection(selection, atoms, molecules, entities);
  }

  std::sort(entities.begin(), entities.end(), [](const auto &a, const auto &b) {
    auto get_sort_key = [](const auto &entity) {
      return std::visit(
          [](const auto &e) {
            return std::make_pair(typeid(e).hash_code(), e.index);
          },
          entity);
    };
    return get_sort_key(a) < get_sort_key(b);
  });

  entities.erase(
      std::unique(entities.begin(), entities.end(),
                  [](const auto &a, const auto &b) {
                    auto is_same_type_and_index = [](const auto &a_entity,
                                                     const auto &b_entity) {
                      return std::visit(
                          [&](const auto &a_val, const auto &b_val) {
                            using TypeA = std::decay_t<decltype(a_val)>;
                            using TypeB = std::decay_t<decltype(b_val)>;
                            if constexpr (std::is_same_v<TypeA, TypeB>) {
                              return a_val.index == b_val.index;
                            }
                            return false;
                          },
                          a_entity, b_entity);
                    };
                    return is_same_type_and_index(a, b);
                  }),
      entities.end());
}

void update_entities_with_positions(std::vector<core::EntityVariant> &entities,
                                    const std::vector<Atom> &atoms,
                                    const std::vector<Molecule> &molecules) {
  for (core::EntityVariant &entity : entities) {
    std::visit(
        [&](auto &e) {
          using ent_type = std::decay_t<decltype(e)>;
          if constexpr (std::is_same_v<ent_type, core::Atom>) {
            e = atoms[e.index];
          } else if constexpr (std::is_same_v<ent_type, core::Molecule>) {
            e = molecules[e.index];
          }
        },
        entity);
  }
}

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
