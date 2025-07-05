#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <fstream>
#include <trajan/core/trajectory.h>
#include <trajan/io/selection.h>

namespace fs = std::filesystem;

using trajan::core::Trajectory;
using trajan::io::AtomIndexSelection;
using trajan::io::AtomTypeSelection;
using trajan::io::MoleculeIndexSelection;
using trajan::io::MoleculeTypeSelection;
using trajan::io::SelectionParser;

fs::path CURRENT = __FILE__;
fs::path EXAMPLES_DIR = CURRENT.parent_path().parent_path() / "examples";

TEST_CASE("Selection Parser - Index Selection", "[selection][index]") {
  SECTION("Single index") {
    auto result = SelectionParser::parse("i1");
    REQUIRE(result.has_value());

    const auto *selection = std::get_if<AtomIndexSelection>(&*result);
    REQUIRE(selection != nullptr);
    REQUIRE(selection->data.size() == 1);
    CHECK(selection->data[0] == 1);
  }

  SECTION("Multiple indices") {
    auto result = SelectionParser::parse("i1,2,3");
    REQUIRE(result.has_value());

    const auto *selection = std::get_if<AtomIndexSelection>(&*result);
    REQUIRE(selection != nullptr);
    REQUIRE(selection->data.size() == 3);
    CHECK(selection->data == std::vector{1, 2, 3});
  }

  SECTION("Index range") {
    auto result = SelectionParser::parse("i1-3");
    REQUIRE(result.has_value());

    const auto *selection = std::get_if<AtomIndexSelection>(&*result);
    REQUIRE(selection != nullptr);
    REQUIRE(selection->data.size() == 3);
    CHECK(selection->data == std::vector{1, 2, 3});
  }

  SECTION("Mixed indices and ranges") {
    auto result = SelectionParser::parse("i1,2-4,6");
    REQUIRE(result.has_value());

    const auto *selection = std::get_if<AtomIndexSelection>(&*result);
    REQUIRE(selection != nullptr);
    REQUIRE(selection->data.size() == 5);
    CHECK(selection->data == std::vector{1, 2, 3, 4, 6});
  }
}

TEST_CASE("Selection Parser - Atom Type Selection", "[selection][atom]") {
  SECTION("Single atom type") {
    auto result = SelectionParser::parse("aC");
    REQUIRE(result.has_value());

    const auto *selection = std::get_if<AtomTypeSelection>(&*result);
    REQUIRE(selection != nullptr);
    REQUIRE(selection->data.size() == 1);
    CHECK(selection->data[0] == "C");
  }

  SECTION("Multiple atom types") {
    auto result = SelectionParser::parse("aC,N,O");
    REQUIRE(result.has_value());

    const auto *selection = std::get_if<AtomTypeSelection>(&*result);
    REQUIRE(selection != nullptr);
    REQUIRE(selection->data.size() == 3);
    CHECK(selection->data == std::vector<std::string>{"C", "N", "O"});
  }

  SECTION("Atom types with underscores") {
    auto result = SelectionParser::parse("aCA_1,CB_2");
    REQUIRE(result.has_value());

    const auto *selection = std::get_if<AtomTypeSelection>(&*result);
    REQUIRE(selection != nullptr);
    REQUIRE(selection->data.size() == 2);
    CHECK(selection->data == std::vector<std::string>{"CA_1", "CB_2"});
  }
}

TEST_CASE("Selection Parser - Molecule Index Selection",
          "[selection][molecule]") {
  SECTION("Single molecule") {
    auto result = SelectionParser::parse("j1");
    REQUIRE(result.has_value());

    const auto *selection = std::get_if<MoleculeIndexSelection>(&*result);
    REQUIRE(selection != nullptr);
    REQUIRE(selection->data.size() == 1);
    CHECK(selection->data[0] == 1);
  }

  SECTION("Molecule range") {
    auto result = SelectionParser::parse("j1-3");
    REQUIRE(result.has_value());

    const auto *selection = std::get_if<MoleculeIndexSelection>(&*result);
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
         {"aC-N", "Range not allowed for atom types"},
         {"i-1", "Negative index"},
         {"j-1", "Negative molecule ID"},
         {"aC@#", "Invalid characters in atom type"},
         {"", "Empty input"}}));

    INFO("Testing: " << description);
    CHECK_FALSE(SelectionParser::parse(input).has_value());
  }
}

TEST_CASE("Selection Parser - Edge Cases", "[selection][edge]") {
  SECTION("Whitespace handling") {
    auto result = SelectionParser::parse("i1, 2,  3");
    REQUIRE(result.has_value());

    const auto *selection = std::get_if<AtomIndexSelection>(&*result);
    REQUIRE(selection != nullptr);
    REQUIRE(selection->data.size() == 3);
    CHECK(selection->data == std::vector{1, 2, 3});
  }

  SECTION("Duplicate values are removed") {
    auto result = SelectionParser::parse("i1,2,2,3,1");
    REQUIRE(result.has_value());

    const auto *selection = std::get_if<AtomIndexSelection>(&*result);
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
                 TestCase{"j1-3,2-4", 4},   // Overlapping molecule ranges
                 TestCase{"aC,C,O,O,N", 3}  // Duplicate atom types
        );

    INFO("Testing input: " << test_case.input);
    auto result = SelectionParser::parse(test_case.input);
    REQUIRE(result.has_value());

    std::visit(
        [&test_case](const auto &selection) {
          using T = std::decay_t<decltype(selection)>;
          if constexpr (std::is_same_v<T, AtomIndexSelection>) {
            CHECK(selection.data.size() == test_case.expected_size);
          } else if constexpr (std::is_same_v<T, AtomTypeSelection>) {
            CHECK(selection.data.size() == test_case.expected_size);
          } else if constexpr (std::is_same_v<T, MoleculeIndexSelection>) {
            CHECK(selection.data.size() == test_case.expected_size);
          }
        },
        *result);
  }
}

TEST_CASE("PDB Read/Write", "[io][pdb]") {
  fs::path temp_pdb = EXAMPLES_DIR / "test_write.pdb";

  Trajectory traj_read;
  std::vector<fs::path> files = {EXAMPLES_DIR / "coord.pdb"};
  traj_read.load_files(files);
  REQUIRE(traj_read.next_frame());

  traj_read.set_output_file(temp_pdb);
  traj_read.write_frame();

  Trajectory traj_write;
  std::vector<fs::path> written_files = {temp_pdb};
  traj_write.load_files(written_files);
  REQUIRE(traj_write.next_frame());

  std::ifstream temp_file(temp_pdb);
  std::string temp_content((std::istreambuf_iterator<char>(temp_file)),
                           std::istreambuf_iterator<char>());
  INFO("Temporary PDB content:\n" << temp_content);
  REQUIRE(!temp_content.empty());

  Trajectory traj_read_original;
  std::vector<fs::path> original_files = {EXAMPLES_DIR / "coord.pdb"};
  traj_read_original.load_files(original_files);
  REQUIRE(traj_read_original.next_frame());

  const auto &atoms_read = traj_read_original.atoms();
  const auto &atoms_write = traj_write.atoms();
  REQUIRE(atoms_read.size() == atoms_write.size());

  for (size_t i = 0; i < atoms_read.size(); ++i) {
    REQUIRE_THAT(atoms_read[i].x,
                 Catch::Matchers::WithinAbs(atoms_write[i].x, 1e-3));
    REQUIRE_THAT(atoms_read[i].y,
                 Catch::Matchers::WithinAbs(atoms_write[i].y, 1e-3));
    REQUIRE_THAT(atoms_read[i].z,
                 Catch::Matchers::WithinAbs(atoms_write[i].z, 1e-3));
  }

  const auto &uc_read = traj_read_original.unit_cell();
  const auto &uc_write = traj_write.unit_cell();
  REQUIRE_THAT(uc_read.a(), Catch::Matchers::WithinAbs(uc_write.a(), 1e-3));
  REQUIRE_THAT(uc_read.b(), Catch::Matchers::WithinAbs(uc_write.b(), 1e-3));
  REQUIRE_THAT(uc_read.c(), Catch::Matchers::WithinAbs(uc_write.c(), 1e-3));
  REQUIRE_THAT(uc_read.alpha(),
               Catch::Matchers::WithinAbs(uc_write.alpha(), 1e-3));
  REQUIRE_THAT(uc_read.beta(),
               Catch::Matchers::WithinAbs(uc_write.beta(), 1e-3));
  REQUIRE_THAT(uc_read.gamma(),
               Catch::Matchers::WithinAbs(uc_write.gamma(), 1e-3));

  fs::remove(temp_pdb);
}

TEST_CASE("DCD Read/Write", "[io][dcd]") {
  fs::path temp_dcd = EXAMPLES_DIR / "test_write.dcd";

  Trajectory traj_read;
  std::vector<fs::path> files = {EXAMPLES_DIR / "coord.pdb",
                                 EXAMPLES_DIR / "traj.dcd"};
  traj_read.load_files(files);

  traj_read.set_output_file(temp_dcd);

  while (traj_read.next_frame()) {
    if (traj_read.current_frame_index() == 1) {
      continue;
    }
    traj_read.write_frame();
  }

  Trajectory traj_write;
  std::vector<fs::path> written_files = {EXAMPLES_DIR / "coord.pdb", temp_dcd};
  traj_write.load_files(written_files);

  std::vector<std::vector<trajan::core::Atom>> all_atoms_write;
  while (traj_write.next_frame()) {
    all_atoms_write.push_back(traj_write.atoms());
  }

  Trajectory traj_read_original;
  std::vector<fs::path> original_files = {EXAMPLES_DIR / "coord.pdb",
                                          EXAMPLES_DIR / "traj.dcd"};
  traj_read_original.load_files(original_files);

  std::vector<std::vector<trajan::core::Atom>> all_atoms_read;
  while (traj_read_original.next_frame()) {
    all_atoms_read.push_back(traj_read_original.atoms());
  }

  REQUIRE(all_atoms_read.size() == all_atoms_write.size());
  double tol = 1e-6;
  for (size_t i = 0; i < all_atoms_read.size(); ++i) {
    const auto &atoms_read = all_atoms_read[i];
    const auto &atoms_write = all_atoms_write[i];
    REQUIRE(atoms_read.size() == atoms_write.size());
    for (size_t j = 0; j < atoms_read.size(); ++j) {
      REQUIRE_THAT(atoms_read[j].x,
                   Catch::Matchers::WithinAbs(atoms_write[j].x, tol));
      REQUIRE_THAT(atoms_read[j].y,
                   Catch::Matchers::WithinAbs(atoms_write[j].y, tol));
      REQUIRE_THAT(atoms_read[j].z,
                   Catch::Matchers::WithinAbs(atoms_write[j].z, tol));
    }
  }

  fs::remove(temp_dcd);
}
