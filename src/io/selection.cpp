#include "trajan/core/neigh.h"
#include <algorithm>
#include <ankerl/unordered_dense.h>
#include <cctype>
#include <occ/core/molecule.h>
#include <stdexcept>
#include <string>
#include <trajan/core/util.h>
#include <trajan/io/selection.h>
#include <unordered_set>

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

  // Flush the accumulated token into results. Returns false on parse failure.
  auto flush_token = [&]() -> bool {
    if (current_token.empty() || current_prefix == '\0')
      return true;
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
    default:
      return false;
    }
    if (!result)
      return false;
    results.push_back(*result);
    current_token.clear();
    return true;
  };

  for (size_t i = 0; i < input.length(); i++) {
    char c = input[i];

    // A prefix character is only recognised at position 0 or after a comma.
    if ((c == SelectionTraits<AtomIndexSelection>::prefix ||
         c == SelectionTraits<AtomTypeSelection>::prefix ||
         c == SelectionTraits<MoleculeIndexSelection>::prefix ||
         c == SelectionTraits<MoleculeTypeSelection>::prefix) &&
        (i == 0 || input[i - 1] == ',')) {

      if (!flush_token())
        return std::nullopt;
      current_prefix = c;
    } else {
      current_token += c;
    }
  }

  if (!flush_token())
    return std::nullopt;

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
  const size_t entities_before = entities.size();
  std::visit(
      [&](const auto &sel) {
        using SelType = std::decay_t<decltype(sel)>;
        trajan::log::debug("Processing selection '{}'", sel.raw_input());
        trajan::log::debug("Selection identified as '{}'", sel.name());

        if constexpr (std::is_same_v<SelType, io::AtomIndexSelection>) {
          trajan::log::debug(
              "The following atoms have been idenitified from {}:", sel.name());
          trajan::log::debug(" NB: atoms are only printed if log level is 4.");
          // sel.data is sorted; use binary search for O(n log k) vs O(nk)
          for (const Atom &atom : atoms) {
            if (std::binary_search(sel.data.begin(), sel.data.end(),
                                   atom.index)) {
              entities.push_back(atom);
              trajan::log::trace(" {}", atom.repr());
            }
          }
        } else if constexpr (std::is_same_v<SelType, io::AtomTypeSelection>) {
          trajan::log::debug(
              "The following atoms have been idenitified from {}:", sel.name());
          trajan::log::debug(" NB: atoms are only printed if log level is 4.");
          // Build a set once for O(n) lookup across all atoms
          const std::unordered_set<std::string> type_set(sel.data.begin(),
                                                         sel.data.end());
          for (const Atom &atom : atoms) {
            if (type_set.count(atom.type)) {
              entities.push_back(atom);
              trajan::log::trace(" {}", atom.repr());
            }
          }
        } else if constexpr (std::is_same_v<SelType,
                                            io::MoleculeIndexSelection>) {
          trajan::log::debug(
              "The following molecules have been idenitified from {}:",
              sel.name());
          trajan::log::debug(" NB: atoms are only printed if log level is 4.");
          // sel.data is sorted; use binary search
          for (const Molecule &molecule : molecules) {
            if (std::binary_search(sel.data.begin(), sel.data.end(),
                                   molecule.index)) {
              entities.push_back(molecule);
              trajan::log::trace("  {}", molecule.repr());
            }
          }
        } else if constexpr (std::is_same_v<SelType,
                                            io::MoleculeTypeSelection>) {
          trajan::log::debug(
              "The following molecules have been idenitified from {}:",
              sel.name());
          trajan::log::debug(" NB: atoms are only printed if log level is 4.");
          // Build a set once for O(n) lookup
          const std::unordered_set<std::string> type_set(sel.data.begin(),
                                                         sel.data.end());
          for (const Molecule &molecule : molecules) {
            if (type_set.count(molecule.type)) {
              entities.push_back(molecule);
              trajan::log::trace("  {}", molecule.repr());
            }
          }
        }

        // Check that THIS criterion matched something (not just that the
        // accumulated entities vector is non-empty from a previous criterion).
        if (entities.size() == entities_before) {
          throw std::runtime_error(fmt::format(
              "No entities found for selection '{}'", sel.raw_input()));
        }
        trajan::log::debug("Identified {} entities in selection '{}'",
                           entities.size() - entities_before, sel.raw_input());
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
    throw std::invalid_argument(fmt::format(
        "Invalid selection: '{}'. "
        "Expected format: i<int>[,<int>|-<int>]... | a<type>[,<type>]... | "
        "j<int>[,<int>|-<int>]... | m<type>[,<type>]...",
        input));
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


// ── Tokenizer ─────────────────────────────────────────────────────────────────

enum class TokKind {
  LParen,
  RParen,
  And,
  Or,
  Not,
  Selector, // e.g. "aO,N", "i1,2-5"
  Coord,    // x, y, or z (as position axis)
  Comp,     // <, <=, >, >=
  Number,
  End
};

struct Token {
  TokKind kind;
  std::string value;
};

static std::vector<Token> tokenize(const std::string &input) {
  std::vector<Token> tokens;
  size_t i = 0;
  const size_t n = input.size();

  auto skip_ws = [&]() {
    while (i < n && std::isspace(static_cast<unsigned char>(input[i])))
      ++i;
  };

  // Case-insensitive prefix match at position i.
  auto match_kw = [&](const char *kw) -> bool {
    size_t len = std::strlen(kw);
    if (i + len > n)
      return false;
    for (size_t k = 0; k < len; ++k) {
      if (std::tolower(static_cast<unsigned char>(input[i + k])) != kw[k])
        return false;
    }
    // Must be followed by whitespace, '(', or end of string.
    size_t after = i + len;
    if (after < n && !std::isspace(static_cast<unsigned char>(input[after])) &&
        input[after] != '(')
      return false;
    return true;
  };

  while (i < n) {
    skip_ws();
    if (i >= n)
      break;

    char c = input[i];

    if (c == '(') {
      tokens.push_back({TokKind::LParen, "("});
      ++i;
      continue;
    }
    if (c == ')') {
      tokens.push_back({TokKind::RParen, ")"});
      ++i;
      continue;
    }

    // Comparators (check two-char before one-char)
    if (c == '>' || c == '<') {
      if (i + 1 < n && input[i + 1] == '=') {
        tokens.push_back({TokKind::Comp, std::string(1, c) + "="});
        i += 2;
      } else {
        tokens.push_back({TokKind::Comp, std::string(1, c)});
        ++i;
      }
      continue;
    }

    // Keywords: and, or, not
    if (match_kw("and")) {
      tokens.push_back({TokKind::And, "and"});
      i += 3;
      continue;
    }
    if (match_kw("or")) {
      tokens.push_back({TokKind::Or, "or"});
      i += 2;
      continue;
    }
    if (match_kw("not")) {
      tokens.push_back({TokKind::Not, "not"});
      i += 3;
      continue;
    }

    // Selector: starts with i, a, j, m followed by non-whitespace non-paren
    if (c == 'i' || c == 'a' || c == 'j' || c == 'm') {
      // Peek ahead: if followed by alphanumeric or comma or hyphen it's a
      // selector. A bare 'a' before a space (e.g. in "a and b") would be
      // ambiguous — treat as selector only if followed by a non-whitespace
      // non-paren character.
      size_t j = i + 1;
      if (j < n && !std::isspace(static_cast<unsigned char>(input[j])) &&
          input[j] != '(' && input[j] != ')') {
        // Consume until whitespace, paren, or comparator
        while (j < n && !std::isspace(static_cast<unsigned char>(input[j])) &&
               input[j] != '(' && input[j] != ')' && input[j] != '<' &&
               input[j] != '>')
          ++j;
        tokens.push_back({TokKind::Selector, input.substr(i, j - i)});
        i = j;
        continue;
      }
    }

    // Coordinate axis: standalone x, y, z
    if ((c == 'x' || c == 'y' || c == 'z') &&
        (i + 1 >= n || std::isspace(static_cast<unsigned char>(input[i + 1])) ||
         input[i + 1] == '<' || input[i + 1] == '>')) {
      tokens.push_back({TokKind::Coord, std::string(1, c)});
      ++i;
      continue;
    }

    // Number (possibly negative)
    if (std::isdigit(static_cast<unsigned char>(c)) ||
        (c == '-' && i + 1 < n &&
         std::isdigit(static_cast<unsigned char>(input[i + 1])))) {
      size_t j = i;
      if (input[j] == '-')
        ++j;
      while (j < n && (std::isdigit(static_cast<unsigned char>(input[j])) ||
                       input[j] == '.'))
        ++j;
      tokens.push_back({TokKind::Number, input.substr(i, j - i)});
      i = j;
      continue;
    }

    throw std::invalid_argument(
        fmt::format("Unexpected character '{}' in selection at position {}", c,
                    i));
  }

  tokens.push_back({TokKind::End, ""});
  return tokens;
}

// ── Recursive-descent parser ──────────────────────────────────────────────────

struct ParseState {
  const std::vector<Token> &tokens;
  size_t pos{0};

  const Token &peek() const { return tokens[pos]; }
  Token consume() { return tokens[pos++]; }
  bool at(TokKind k) const { return tokens[pos].kind == k; }
};

static SelectionExpr parse_expr_impl(ParseState &ps);
static SelectionExpr parse_or(ParseState &ps);
static SelectionExpr parse_and(ParseState &ps);
static SelectionExpr parse_not(ParseState &ps);
static SelectionExpr parse_primary(ParseState &ps);

static SelectionExpr parse_expr_impl(ParseState &ps) { return parse_or(ps); }

static SelectionExpr parse_or(ParseState &ps) {
  SelectionExpr lhs = parse_and(ps);
  while (ps.at(TokKind::Or)) {
    ps.consume();
    SelectionExpr rhs = parse_and(ps);
    OrExpr node;
    node.lhs = std::make_shared<SelectionExpr>(std::move(lhs));
    node.rhs = std::make_shared<SelectionExpr>(std::move(rhs));
    lhs = SelectionExpr{std::move(node)};
  }
  return lhs;
}

static SelectionExpr parse_and(ParseState &ps) {
  SelectionExpr lhs = parse_not(ps);
  while (ps.at(TokKind::And)) {
    ps.consume();
    SelectionExpr rhs = parse_not(ps);
    AndExpr node;
    node.lhs = std::make_shared<SelectionExpr>(std::move(lhs));
    node.rhs = std::make_shared<SelectionExpr>(std::move(rhs));
    lhs = SelectionExpr{std::move(node)};
  }
  return lhs;
}

static SelectionExpr parse_not(ParseState &ps) {
  if (ps.at(TokKind::Not)) {
    ps.consume();
    SelectionExpr operand = parse_primary(ps);
    NotExpr node;
    node.operand = std::make_shared<SelectionExpr>(std::move(operand));
    return SelectionExpr{std::move(node)};
  }
  return parse_primary(ps);
}

static SelectionExpr parse_primary(ParseState &ps) {
  if (ps.at(TokKind::LParen)) {
    ps.consume();
    SelectionExpr inner = parse_expr_impl(ps);
    if (!ps.at(TokKind::RParen)) {
      throw std::invalid_argument("Expected ')' in selection expression");
    }
    ps.consume();
    return inner;
  }

  if (ps.at(TokKind::Selector)) {
    Token tok = ps.consume();
    auto result = SelectionParser::parse(tok.value);
    if (!result) {
      throw std::invalid_argument(
          fmt::format("Invalid selector token '{}'", tok.value));
    }
    return SelectionExpr{CriteriaLeaf{std::move(*result)}};
  }

  if (ps.at(TokKind::Coord)) {
    Token coord_tok = ps.consume();
    if (!ps.at(TokKind::Comp)) {
      throw std::invalid_argument(
          fmt::format("Expected comparator after coordinate '{}'",
                      coord_tok.value));
    }
    Token comp_tok = ps.consume();
    if (!ps.at(TokKind::Number)) {
      throw std::invalid_argument(
          fmt::format("Expected number after comparator '{}'", comp_tok.value));
    }
    Token num_tok = ps.consume();

    Axis axis;
    if (coord_tok.value == "x")
      axis = Axis::X;
    else if (coord_tok.value == "y")
      axis = Axis::Y;
    else
      axis = Axis::Z;

    CompOp op;
    if (comp_tok.value == "<")
      op = CompOp::LT;
    else if (comp_tok.value == "<=")
      op = CompOp::LE;
    else if (comp_tok.value == ">")
      op = CompOp::GT;
    else
      op = CompOp::GE;

    double val = std::stod(num_tok.value);
    return SelectionExpr{PositionLeaf{{axis, op, val}}};
  }

  throw std::invalid_argument(
      fmt::format("Unexpected token '{}' in selection expression",
                  ps.peek().value));
}

// Check whether an input string uses any extended-expression syntax.
// If not, we can delegate to the simpler SelectionParser::parse() path.
static bool is_simple_selection(const std::string &input) {
  for (size_t i = 0; i < input.size(); ++i) {
    char c = input[i];
    if (c == '(' || c == ')' || c == '<' || c == '>')
      return false;
    // Check for keywords at word boundaries
    auto at_word = [&](const char *kw) {
      size_t len = std::strlen(kw);
      if (i + len > input.size())
        return false;
      for (size_t k = 0; k < len; ++k) {
        if (std::tolower(static_cast<unsigned char>(input[i + k])) != kw[k])
          return false;
      }
      size_t after = i + len;
      bool before_ok = (i == 0 || std::isspace(static_cast<unsigned char>(input[i - 1])));
      bool after_ok = (after >= input.size() || std::isspace(static_cast<unsigned char>(input[after])));
      return before_ok && after_ok;
    };
    if (at_word("and") || at_word("or") || at_word("not"))
      return false;
  }
  return true;
}

SelectionExpr SelectionParser::parse_expr(const std::string &input) {
  if (input.empty()) {
    throw std::invalid_argument("Empty selection expression");
  }

  // Fast path: no extended syntax → wrap existing parse() result.
  if (is_simple_selection(input)) {
    auto result = SelectionParser::parse(input);
    if (!result) {
      throw std::invalid_argument(
          fmt::format("Invalid selection: '{}'", input));
    }
    return SelectionExpr{CriteriaLeaf{std::move(*result)}};
  }

  auto tokens = tokenize(input);
  ParseState ps{tokens};
  SelectionExpr expr = parse_expr_impl(ps);
  if (!ps.at(TokKind::End)) {
    throw std::invalid_argument(
        fmt::format("Unexpected trailing token '{}' in selection expression",
                    ps.peek().value));
  }
  return expr;
}

// ── Evaluator ─────────────────────────────────────────────────────────────────

struct EvalResult {
  enum class Kind { Atom, Molecule, PositionOnly } kind;
  ankerl::unordered_dense::set<size_t> indices;
  // Only populated when kind == PositionOnly; all predicates must hold (AND).
  std::vector<PositionPredicate> position_preds;
};

static occ::core::Molecule::Origin to_occ_origin(MolOrigin o) {
  return o == MolOrigin::CenterOfMass ? occ::core::Molecule::CenterOfMass
                                      : occ::core::Molecule::Centroid;
}

static double get_coord(const Atom &atom, Axis axis) {
  switch (axis) {
  case Axis::X: return atom.x;
  case Axis::Y: return atom.y;
  case Axis::Z: return atom.z;
  }
  return 0.0;
}

static double get_coord(const occ::Vec3 &pos, Axis axis) {
  switch (axis) {
  case Axis::X: return pos[0];
  case Axis::Y: return pos[1];
  case Axis::Z: return pos[2];
  }
  return 0.0;
}

// Apply a list of position predicates (all must hold) to a typed entity set.
static EvalResult apply_position_filter(const std::vector<PositionPredicate> &preds,
                                         const EvalResult &entities,
                                         const std::vector<Atom> &atoms,
                                         const std::vector<Molecule> &molecules,
                                         MolOrigin mol_origin) {
  EvalResult result;
  result.kind = entities.kind;
  for (size_t idx : entities.indices) {
    bool pass = true;
    if (entities.kind == EvalResult::Kind::Atom) {
      for (const auto &p : preds) {
        if (!p.evaluate(get_coord(atoms[idx], p.axis))) {
          pass = false;
          break;
        }
      }
    } else {
      occ::Vec3 pos = molecules[idx].position(to_occ_origin(mol_origin));
      for (const auto &p : preds) {
        if (!p.evaluate(get_coord(pos, p.axis))) {
          pass = false;
          break;
        }
      }
    }
    if (pass)
      result.indices.insert(idx);
  }
  return result;
}

// Evaluate PositionOnly as all atoms satisfying all predicates.
static EvalResult resolve_position_only(const EvalResult &pos_result,
                                         const std::vector<Atom> &atoms) {
  EvalResult result;
  result.kind = EvalResult::Kind::Atom;
  for (size_t i = 0; i < atoms.size(); ++i) {
    bool pass = true;
    for (const auto &p : pos_result.position_preds) {
      if (!p.evaluate(get_coord(atoms[i], p.axis))) {
        pass = false;
        break;
      }
    }
    if (pass)
      result.indices.insert(i);
  }
  return result;
}

static EvalResult eval_node(const SelectionExpr &expr,
                             const std::vector<Atom> &atoms,
                             const std::vector<Molecule> &molecules,
                             MolOrigin mol_origin);

static EvalResult eval_node(const SelectionExpr &expr,
                             const std::vector<Atom> &atoms,
                             const std::vector<Molecule> &molecules,
                             MolOrigin mol_origin) {
  return std::visit(
      [&](const auto &node) -> EvalResult {
        using T = std::decay_t<decltype(node)>;

        if constexpr (std::is_same_v<T, CriteriaLeaf>) {
          std::vector<core::EntityVariant> ents;
          process_selection(node.criteria, atoms, molecules, ents);
          EvalResult r;
          // Determine kind from criteria types.
          bool has_mol = false;
          for (const auto &c : node.criteria) {
            if (std::holds_alternative<MoleculeIndexSelection>(c) ||
                std::holds_alternative<MoleculeTypeSelection>(c)) {
              has_mol = true;
              break;
            }
          }
          r.kind = has_mol ? EvalResult::Kind::Molecule : EvalResult::Kind::Atom;
          for (const auto &ev : ents) {
            std::visit(
                [&](const auto &e) { r.indices.insert(static_cast<size_t>(e.index)); },
                ev);
          }
          return r;

        } else if constexpr (std::is_same_v<T, PositionLeaf>) {
          // Return a PositionOnly result — evaluated in context by the parent.
          EvalResult r;
          r.kind = EvalResult::Kind::PositionOnly;
          r.position_preds.push_back(node.pred);
          return r;

        } else if constexpr (std::is_same_v<T, AndExpr>) {
          EvalResult lhs = eval_node(*node.lhs, atoms, molecules, mol_origin);
          EvalResult rhs = eval_node(*node.rhs, atoms, molecules, mol_origin);

          bool lhs_pos = lhs.kind == EvalResult::Kind::PositionOnly;
          bool rhs_pos = rhs.kind == EvalResult::Kind::PositionOnly;

          if (lhs_pos && rhs_pos) {
            // Combine predicates — still PositionOnly.
            EvalResult r;
            r.kind = EvalResult::Kind::PositionOnly;
            r.position_preds.insert(r.position_preds.end(),
                                    lhs.position_preds.begin(),
                                    lhs.position_preds.end());
            r.position_preds.insert(r.position_preds.end(),
                                    rhs.position_preds.begin(),
                                    rhs.position_preds.end());
            return r;
          }
          if (lhs_pos && !rhs_pos) {
            return apply_position_filter(lhs.position_preds, rhs, atoms,
                                          molecules, mol_origin);
          }
          if (!lhs_pos && rhs_pos) {
            return apply_position_filter(rhs.position_preds, lhs, atoms,
                                          molecules, mol_origin);
          }
          // Both typed — must be same kind.
          if (lhs.kind != rhs.kind) {
            throw std::invalid_argument(
                "Cannot AND atom selection with molecule selection");
          }
          EvalResult r;
          r.kind = lhs.kind;
          for (size_t idx : lhs.indices) {
            if (rhs.indices.count(idx))
              r.indices.insert(idx);
          }
          return r;

        } else if constexpr (std::is_same_v<T, OrExpr>) {
          EvalResult lhs = eval_node(*node.lhs, atoms, molecules, mol_origin);
          EvalResult rhs = eval_node(*node.rhs, atoms, molecules, mol_origin);

          // Resolve any PositionOnly to Atom context before OR.
          if (lhs.kind == EvalResult::Kind::PositionOnly)
            lhs = resolve_position_only(lhs, atoms);
          if (rhs.kind == EvalResult::Kind::PositionOnly)
            rhs = resolve_position_only(rhs, atoms);

          if (lhs.kind != rhs.kind) {
            throw std::invalid_argument(
                "Cannot OR atom selection with molecule selection");
          }
          EvalResult r;
          r.kind = lhs.kind;
          r.indices = lhs.indices;
          for (size_t idx : rhs.indices)
            r.indices.insert(idx);
          return r;

        } else if constexpr (std::is_same_v<T, NotExpr>) {
          EvalResult operand =
              eval_node(*node.operand, atoms, molecules, mol_origin);
          if (operand.kind == EvalResult::Kind::PositionOnly)
            operand = resolve_position_only(operand, atoms);

          EvalResult r;
          r.kind = operand.kind;
          size_t total = (operand.kind == EvalResult::Kind::Atom)
                             ? atoms.size()
                             : molecules.size();
          for (size_t i = 0; i < total; ++i) {
            if (!operand.indices.count(i))
              r.indices.insert(i);
          }
          return r;
        }

        // Unreachable.
        return EvalResult{};
      },
      expr.node);
}

void process_selection(const SelectionExpr &expr,
                       const std::vector<Atom> &atoms,
                       const std::vector<Molecule> &molecules,
                       std::vector<core::EntityVariant> &entities,
                       MolOrigin mol_origin) {
  EvalResult result = eval_node(expr, atoms, molecules, mol_origin);

  // Resolve any top-level PositionOnly to Atom context.
  if (result.kind == EvalResult::Kind::PositionOnly)
    result = resolve_position_only(result, atoms);

  if (result.kind == EvalResult::Kind::Atom) {
    for (size_t idx : result.indices)
      entities.push_back(atoms[idx]);
  } else {
    for (size_t idx : result.indices)
      entities.push_back(molecules[idx]);
  }
}

SelectionExpr selection_expr_validator(const std::string &input) {
  trajan::log::debug("Raw selection expression: '{}'", input);
  return SelectionParser::parse_expr(input);
}

} // namespace trajan::io
