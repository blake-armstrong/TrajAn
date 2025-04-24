#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <trajan/io/selection.h>

using trajan::io::AtomTypeSelection;
using trajan::io::IndexSelection;
using trajan::io::MoleculeSelection;
using trajan::io::SelectionParser;

TEST_CASE("Selection Parser - Index Selection", "[selection][index]") {
  SECTION("Single index") {
    auto result = SelectionParser::parse("i1");
    REQUIRE(result.has_value());

    const auto *selection = std::get_if<IndexSelection>(&*result);
    REQUIRE(selection != nullptr);
    REQUIRE(selection->data.size() == 1);
    CHECK(selection->data[0] == 1);
  }

  SECTION("Multiple indices") {
    auto result = SelectionParser::parse("i1,2,3");
    REQUIRE(result.has_value());

    const auto *selection = std::get_if<IndexSelection>(&*result);
    REQUIRE(selection != nullptr);
    REQUIRE(selection->data.size() == 3);
    CHECK(selection->data == std::vector{1, 2, 3});
  }

  SECTION("Index range") {
    auto result = SelectionParser::parse("i1-3");
    REQUIRE(result.has_value());

    const auto *selection = std::get_if<IndexSelection>(&*result);
    REQUIRE(selection != nullptr);
    REQUIRE(selection->data.size() == 3);
    CHECK(selection->data == std::vector{1, 2, 3});
  }

  SECTION("Mixed indices and ranges") {
    auto result = SelectionParser::parse("i1,2-4,6");
    REQUIRE(result.has_value());

    const auto *selection = std::get_if<IndexSelection>(&*result);
    REQUIRE(selection != nullptr);
    REQUIRE(selection->data.size() == 5);
    CHECK(selection->data == std::vector{1, 2, 3, 4, 6});
  }
}

TEST_CASE("Selection Parser - Atom Type Selection", "[selection][atom]") {
  SECTION("Single atom type") {
    auto result = SelectionParser::parse("tC");
    REQUIRE(result.has_value());

    const auto *selection = std::get_if<AtomTypeSelection>(&*result);
    REQUIRE(selection != nullptr);
    REQUIRE(selection->data.size() == 1);
    CHECK(selection->data[0] == "C");
  }

  SECTION("Multiple atom types") {
    auto result = SelectionParser::parse("tC,N,O");
    REQUIRE(result.has_value());

    const auto *selection = std::get_if<AtomTypeSelection>(&*result);
    REQUIRE(selection != nullptr);
    REQUIRE(selection->data.size() == 3);
    CHECK(selection->data == std::vector<std::string>{"C", "N", "O"});
  }

  SECTION("Atom types with underscores") {
    auto result = SelectionParser::parse("tCA_1,CB_2");
    REQUIRE(result.has_value());

    const auto *selection = std::get_if<AtomTypeSelection>(&*result);
    REQUIRE(selection != nullptr);
    REQUIRE(selection->data.size() == 2);
    CHECK(selection->data == std::vector<std::string>{"CA_1", "CB_2"});
  }
}

TEST_CASE("Selection Parser - Molecule Selection", "[selection][molecule]") {
  SECTION("Single molecule") {
    auto result = SelectionParser::parse("m1");
    REQUIRE(result.has_value());

    const auto *selection = std::get_if<MoleculeSelection>(&*result);
    REQUIRE(selection != nullptr);
    REQUIRE(selection->data.size() == 1);
    CHECK(selection->data[0] == 1);
  }

  SECTION("Molecule range") {
    auto result = SelectionParser::parse("m1-3");
    REQUIRE(result.has_value());

    const auto *selection = std::get_if<MoleculeSelection>(&*result);
    REQUIRE(selection != nullptr);
    REQUIRE(selection->data.size() == 3);
    CHECK(selection->data == std::vector{1, 2, 3});
  }
}

TEST_CASE("Selection Parser - Invalid Inputs", "[selection][invalid]") {
  SECTION("Invalid selections should return nullopt") {
    auto [input, description] = GENERATE(table<std::string, std::string>(
        {{"x1,2,3", "Invalid prefix"},
         {"i3-1", "End less than start"},
         {"tC-N", "Range not allowed for atom types"},
         {"i-1", "Negative index"},
         {"m-1", "Negative molecule ID"},
         {"tC@#", "Invalid characters in atom type"},
         {"", "Empty input"}}));

    INFO("Testing: " << description);
    CHECK_FALSE(SelectionParser::parse(input).has_value());
  }
}

TEST_CASE("Selection Parser - Edge Cases", "[selection][edge]") {
  SECTION("Whitespace handling") {
    auto result = SelectionParser::parse("i1, 2,  3");
    REQUIRE(result.has_value());

    const auto *selection = std::get_if<IndexSelection>(&*result);
    REQUIRE(selection != nullptr);
    REQUIRE(selection->data.size() == 3);
    CHECK(selection->data == std::vector{1, 2, 3});
  }

  SECTION("Duplicate values are removed") {
    auto result = SelectionParser::parse("i1,2,2,3,1");
    REQUIRE(result.has_value());

    const auto *selection = std::get_if<IndexSelection>(&*result);
    REQUIRE(selection != nullptr);
    REQUIRE(selection->data.size() == 3);
    CHECK(selection->data == std::vector{1, 2, 3});
  }

  SECTION("Complex mixed cases") {
    struct TestCase {
      std::string input;
      size_t expected_size;
    };

    auto test_case =
        GENERATE(TestCase{"i1-3,2-4,5", 5}, // Overlapping ranges
                 TestCase{"i1,1-3,2,3", 3}, // Duplicates in ranges
                 TestCase{"m1-3,2-4", 4},   // Overlapping molecule ranges
                 TestCase{"tC,C,O,O,N", 3}  // Duplicate atom types
        );

    INFO("Testing input: " << test_case.input);
    auto result = SelectionParser::parse(test_case.input);
    REQUIRE(result.has_value());

    std::visit(
        [&test_case](const auto &selection) {
          using T = std::decay_t<decltype(selection)>;
          if constexpr (std::is_same_v<T, IndexSelection>) {
            CHECK(selection.data.size() == test_case.expected_size);
          } else if constexpr (std::is_same_v<T, AtomTypeSelection>) {
            CHECK(selection.data.size() == test_case.expected_size);
          } else if constexpr (std::is_same_v<T, MoleculeSelection>) {
            CHECK(selection.data.size() == test_case.expected_size);
          }
        },
        *result);
  }
}
