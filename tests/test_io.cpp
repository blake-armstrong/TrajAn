#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <filesystem>
#include <fstream>
#include <occ/core/units.h>
#include <occ/crystal/unitcell.h>
#include <trajan/core/atomgraph.h>
#include <trajan/core/log.h>
#include <trajan/core/molecule.h>
#include <trajan/core/neigh.h>
#include <trajan/core/topology.h>
#include <trajan/core/trajectory.h>
#include <trajan/io/pdb.h>
#include <trajan/io/xtc.h>
#include <trajan/io/xyz.h>
#include <trajan/io/selection.h>
#include <vector>

namespace fs = std::filesystem;

using occ::Mat3N;
using occ::Vec3;

using trajan::core::Atom;
using trajan::core::AtomGraph;
using trajan::core::Bond;
using trajan::core::Molecule;
using trajan::core::Trajectory;

using trajan::io::AtomIndexSelection;
using trajan::io::AtomTypeSelection;
using trajan::io::MoleculeIndexSelection;
using trajan::io::MoleculeTypeSelection;
using trajan::io::SelectionCriteria;
using trajan::io::SelectionParser;

static std::string g_test_data_path;

std::string get_test_data_path() { return g_test_data_path; }

bool has_test_data() {
  return !g_test_data_path.empty() && fs::exists(g_test_data_path);
}

class TestFixture {
public:
  fs::path get_test_file(const fs::path &file) {
    this->require_test_data();
    return fs::path(get_test_data_path()) / file;
  }

  void require_test_data() {
    if (!has_test_data()) {
      SKIP("Test data path not provided or doesn't exist. Use --data-path "
           "<path>");
    }
  }
};

int main(int argc, char *argv[]) {
  Catch::Session session;

  auto cli = session.cli() |
             Catch::Clara::Opt(g_test_data_path, "path")["-d"]["--data-path"](
                 "Path to trajectory test data directory");

  session.cli(cli);

  int result = session.applyCommandLine(argc, argv);
  if (result != 0)
    return result;

  if (g_test_data_path.empty()) {
    trajan::log::info(
        "No test data path provided. File-based tests will be skipped.");
    trajan::log::info("Use --data-path <path> to enable all tests.");
  } else if (!fs::exists(g_test_data_path)) {
    trajan::log::info("Warning: Test data path doesn't exist: {}",
                      g_test_data_path);
    trajan::log::info("File-based tests will be skipped.");
    g_test_data_path.clear();
  } else {
    trajan::log::info("Using test data path: {}", g_test_data_path);
  }

  return session.run();
}

// tests that dont require files

TEST_CASE("Selection Parser - Index Selection", "[unit][selection]") {
  SECTION("Single index") {
    auto result = SelectionParser::parse("i1");
    REQUIRE(result.has_value());
    REQUIRE(result.value().size() == 1);

    for (const auto &sel : result.value()) {
      const auto *selection = std::get_if<AtomIndexSelection>(&sel);
      REQUIRE(selection != nullptr);
      REQUIRE(selection->data.size() == 1);
      CHECK(selection->data[0] == 1);
    }
  }

  SECTION("Multiple indices") {
    auto result = SelectionParser::parse("i1,2,3");
    REQUIRE(result.has_value());
    REQUIRE(result.value().size() == 1);

    for (const auto &sel : result.value()) {
      const auto *selection = std::get_if<AtomIndexSelection>(&sel);
      REQUIRE(selection != nullptr);
      REQUIRE(selection->data.size() == 3);
      CHECK(selection->data == std::vector{1, 2, 3});
    }
  }

  SECTION("Index range") {
    auto result = SelectionParser::parse("i1-3");
    REQUIRE(result.has_value());
    REQUIRE(result.value().size() == 1);

    for (const auto &sel : result.value()) {
      const auto *selection = std::get_if<AtomIndexSelection>(&sel);
      REQUIRE(selection != nullptr);
      REQUIRE(selection->data.size() == 3);
      CHECK(selection->data == std::vector{1, 2, 3});
    }
  }

  SECTION("Mixed indices and ranges") {
    auto result = SelectionParser::parse("i1,2-4,6");
    REQUIRE(result.has_value());
    REQUIRE(result.value().size() == 1);

    for (const auto &sel : result.value()) {
      const auto *selection = std::get_if<AtomIndexSelection>(&sel);
      REQUIRE(selection != nullptr);
      REQUIRE(selection->data.size() == 5);
      CHECK(selection->data == std::vector{1, 2, 3, 4, 6});
    }
  }
}

TEST_CASE("Selection Parser - Atom Type Selection", "[unit][selection]") {
  SECTION("Single atom type") {
    auto result = SelectionParser::parse("aC");
    REQUIRE(result.has_value());
    REQUIRE(result.value().size() == 1);

    for (const auto &sel : result.value()) {
      const auto *selection = std::get_if<AtomTypeSelection>(&sel);
      REQUIRE(selection != nullptr);
      REQUIRE(selection->data.size() == 1);
      CHECK(selection->data[0] == "C");
    }
  }

  SECTION("Multiple atom types") {
    auto result = SelectionParser::parse("aC,N,O");
    REQUIRE(result.has_value());
    REQUIRE(result.value().size() == 1);

    for (const auto &sel : result.value()) {
      const auto *selection = std::get_if<AtomTypeSelection>(&sel);
      REQUIRE(selection != nullptr);
      REQUIRE(selection->data.size() == 3);
      CHECK(selection->data == std::vector<std::string>{"C", "N", "O"});
    }
  }

  SECTION("Atom types with underscores") {
    auto result = SelectionParser::parse("aCA_1,CB_2");
    REQUIRE(result.has_value());
    REQUIRE(result.value().size() == 1);

    for (const auto &sel : result.value()) {
      const auto *selection = std::get_if<AtomTypeSelection>(&sel);
      REQUIRE(selection != nullptr);
      REQUIRE(selection->data.size() == 2);
      CHECK(selection->data == std::vector<std::string>{"CA_1", "CB_2"});
    }
  }
}

TEST_CASE("Selection Parser - Molecule Index Selection", "[unit][selection]") {
  SECTION("Single molecule") {
    auto result = SelectionParser::parse("j1");
    REQUIRE(result.has_value());
    REQUIRE(result.value().size() == 1);

    for (const auto &sel : result.value()) {
      const auto *selection = std::get_if<MoleculeIndexSelection>(&sel);
      REQUIRE(selection != nullptr);
      REQUIRE(selection->data.size() == 1);
      CHECK(selection->data[0] == 1);
    }
  }

  SECTION("Molecule range") {
    auto result = SelectionParser::parse("j1-3");
    REQUIRE(result.has_value());
    REQUIRE(result.value().size() == 1);

    for (const auto &sel : result.value()) {
      const auto *selection = std::get_if<MoleculeIndexSelection>(&sel);
      REQUIRE(selection != nullptr);
      REQUIRE(selection->data.size() == 3);
      CHECK(selection->data == std::vector{1, 2, 3});
    }
  }
}

TEST_CASE("Selection Parser - Invalid Inputs", "[unit][selection]") {
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

TEST_CASE("Selection Parser - Edge Cases", "[unit][selection]") {
  SECTION("Whitespace handling") {
    auto result = SelectionParser::parse("i1, 2,  3");
    REQUIRE(result.has_value());
    REQUIRE(result.value().size() == 1);

    for (const auto &sel : result.value()) {
      const auto *selection = std::get_if<AtomIndexSelection>(&sel);
      REQUIRE(selection != nullptr);
      REQUIRE(selection->data.size() == 3);
      CHECK(selection->data == std::vector{1, 2, 3});
    }
  }

  SECTION("Duplicate values are removed") {
    auto result = SelectionParser::parse("i1,2,2,3,1");
    REQUIRE(result.has_value());
    REQUIRE(result.value().size() == 1);

    for (const auto &sel : result.value()) {
      const auto *selection = std::get_if<AtomIndexSelection>(&sel);
      REQUIRE(selection != nullptr);
      REQUIRE(selection->data.size() == 3);
      CHECK(selection->data == std::vector{1, 2, 3});
    }
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
    REQUIRE(result->size() == 1);

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
        result->front());
  }
}

// tests requiring files

TEST_CASE_METHOD(TestFixture, "PDB Read/Write", "[file][io][pdb]") {
  fs::path temp_pdb = "/tmp/test_write.pdb";

  Trajectory traj_read;
  fs::path test_file = this->get_test_file("i00000001.pdb");
  std::vector<fs::path> files = {test_file};
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

  REQUIRE(!temp_content.empty());

  Trajectory traj_read_original;
  traj_read_original.load_files(files);
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

  const auto &uc_read = traj_read_original.unit_cell().value();
  const auto &uc_write = traj_write.unit_cell().value();
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

TEST_CASE_METHOD(TestFixture, "DCD Read/Write", "[file][io][dcd]") {
  fs::path temp_dcd = "/tmp/test_write.dcd";

  Trajectory traj_read;
  fs::path pdb_file = this->get_test_file("i00000001.pdb");
  std::vector<fs::path> files = {pdb_file,
                                 this->get_test_file("i00000002.dcd")};
  traj_read.load_files(files);
  traj_read.set_output_file(temp_dcd);

  while (traj_read.next_frame()) {
    if (traj_read.current_frame_index() == 1) {
      continue;
    }
    traj_read.write_frame();
  }

  Trajectory traj_write;
  std::vector<fs::path> written_files = {pdb_file, temp_dcd};
  traj_write.load_files(written_files);

  std::vector<std::vector<trajan::core::Atom>> all_atoms_write;
  while (traj_write.next_frame()) {
    all_atoms_write.push_back(traj_write.atoms());
  }

  Trajectory traj_read_original;
  traj_read_original.load_files(files);

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

TEST_CASE_METHOD(TestFixture, "PDB Read/Write into memory", "[file][io][pdb]") {

  fs::path temp_pdb = "/tmp/test_write.pdb";

  Trajectory traj_read;
  std::vector<fs::path> files = {this->get_test_file("i00000001.pdb")};
  traj_read.load_files_into_memory(files);
  REQUIRE(traj_read.next_frame());

  traj_read.set_output_file(temp_pdb);
  traj_read.write_frame();

  Trajectory traj_write;
  std::vector<fs::path> written_files = {temp_pdb};
  traj_write.load_files_into_memory(written_files);
  REQUIRE(traj_write.next_frame());

  std::ifstream temp_file(temp_pdb);
  std::string temp_content((std::istreambuf_iterator<char>(temp_file)),
                           std::istreambuf_iterator<char>());
  REQUIRE(!temp_content.empty());

  Trajectory traj_read_original;
  traj_read_original.load_files_into_memory(files);
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

  const auto &uc_read = traj_read_original.unit_cell().value();
  const auto &uc_write = traj_write.unit_cell().value();
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

TEST_CASE_METHOD(TestFixture, "DCD Read/Write into memory", "[file][io][dcd]") {
  fs::path temp_dcd = "/tmp/test_write.dcd";

  Trajectory traj_read;
  fs::path pdb_file = this->get_test_file("i00000001.pdb");
  std::vector<fs::path> files = {pdb_file,
                                 this->get_test_file("i00000002.dcd")};
  traj_read.load_files_into_memory(files);

  traj_read.set_output_file(temp_dcd);

  std::vector<std::vector<trajan::core::Atom>> all_atoms_read;
  while (traj_read.next_frame()) {
    all_atoms_read.push_back(traj_read.atoms());
    if (traj_read.current_frame_index() == 1) {
      continue;
    }
    traj_read.write_frame();
  }

  Trajectory traj_write;
  std::vector<fs::path> written_files = {pdb_file, temp_dcd};
  traj_write.load_files_into_memory(written_files);

  std::vector<std::vector<trajan::core::Atom>> all_atoms_write;
  while (traj_write.next_frame()) {
    all_atoms_write.push_back(traj_write.atoms());
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

// ── XYZ / extended-XYZ tests (self-contained, no external test data needed) ──

static std::vector<trajan::core::Atom>
make_water_atoms() {
  using trajan::core::Atom;
  using trajan::core::Element;
  std::vector<Atom> atoms;
  atoms.emplace_back(occ::Vec3{0.000, 0.000, 0.000}, Element("O", true), 0);
  atoms[0].type = "O";
  atoms.emplace_back(occ::Vec3{0.960, 0.000, 0.000}, Element("H", true), 1);
  atoms[1].type = "H";
  atoms.emplace_back(occ::Vec3{-0.240, 0.927, 0.000}, Element("H", true), 2);
  atoms[2].type = "H";
  return atoms;
}

TEST_CASE("XYZ Write and Read Back - no unit cell", "[unit][io][xyz]") {
  fs::path temp = "/tmp/test_xyz_basic.xyz";

  {
    trajan::io::XYZHandler writer;
    writer.set_file_path(temp);
    writer.initialise(trajan::io::FileHandler::Mode::Write);

    trajan::core::Frame frame;
    frame.set_atoms(make_water_atoms());
    REQUIRE(writer.write_frame(frame));
    writer.finalise();
  }

  trajan::core::Frame frame;
  trajan::io::XYZHandler reader;
  reader.set_file_path(temp);
  reader.initialise(trajan::io::FileHandler::Mode::Read);
  REQUIRE(reader.read_frame(frame));
  reader.finalise();

  REQUIRE(frame.num_atoms() == 3);
  REQUIRE_FALSE(frame.has_unit_cell());

  const auto &atoms = frame.atoms();
  REQUIRE(atoms[0].element.symbol() == "O");
  REQUIRE(atoms[1].element.symbol() == "H");
  REQUIRE(atoms[2].element.symbol() == "H");
  REQUIRE_THAT(atoms[0].x, Catch::Matchers::WithinAbs(0.000, 1e-5));
  REQUIRE_THAT(atoms[1].x, Catch::Matchers::WithinAbs(0.960, 1e-5));
  REQUIRE_THAT(atoms[2].x, Catch::Matchers::WithinAbs(-0.240, 1e-5));
  REQUIRE_THAT(atoms[2].y, Catch::Matchers::WithinAbs(0.927, 1e-5));

  fs::remove(temp);
}

TEST_CASE("XYZ Write and Read Back - extended XYZ with orthogonal Lattice",
          "[unit][io][xyz]") {
  fs::path temp = "/tmp/test_xyz_extended.xyz";
  constexpr double a = 12.3, b = 13.4, c = 14.5;
  const double alpha = occ::units::radians(90.0);
  const double beta  = occ::units::radians(90.0);
  const double gamma = occ::units::radians(90.0);

  {
    trajan::io::XYZHandler writer;
    writer.set_file_path(temp);
    writer.initialise(trajan::io::FileHandler::Mode::Write);

    trajan::core::Frame frame;
    frame.set_atoms(make_water_atoms());
    auto uc = occ::crystal::triclinic_cell(a, b, c, alpha, beta, gamma);
    frame.set_unit_cell(uc);
    REQUIRE(writer.write_frame(frame));
    writer.finalise();
  }

  trajan::core::Frame frame;
  trajan::io::XYZHandler reader;
  reader.set_file_path(temp);
  reader.initialise(trajan::io::FileHandler::Mode::Read);
  REQUIRE(reader.read_frame(frame));
  reader.finalise();

  REQUIRE(frame.num_atoms() == 3);
  REQUIRE(frame.has_unit_cell());

  const auto &uc = frame.unit_cell().value();
  REQUIRE_THAT(uc.a(), Catch::Matchers::WithinAbs(a, 1e-4));
  REQUIRE_THAT(uc.b(), Catch::Matchers::WithinAbs(b, 1e-4));
  REQUIRE_THAT(uc.c(), Catch::Matchers::WithinAbs(c, 1e-4));
  REQUIRE_THAT(occ::units::degrees(uc.alpha()),
               Catch::Matchers::WithinAbs(90.0, 1e-3));
  REQUIRE_THAT(occ::units::degrees(uc.beta()),
               Catch::Matchers::WithinAbs(90.0, 1e-3));
  REQUIRE_THAT(occ::units::degrees(uc.gamma()),
               Catch::Matchers::WithinAbs(90.0, 1e-3));

  fs::remove(temp);
}

TEST_CASE("XYZ Write and Read Back - triclinic cell", "[unit][io][xyz]") {
  fs::path temp = "/tmp/test_xyz_triclinic.xyz";
  constexpr double a = 10.0, b = 11.0, c = 12.0;
  const double alpha = occ::units::radians(80.0);
  const double beta  = occ::units::radians(85.0);
  const double gamma = occ::units::radians(70.0);

  {
    trajan::io::XYZHandler writer;
    writer.set_file_path(temp);
    writer.initialise(trajan::io::FileHandler::Mode::Write);

    trajan::core::Frame frame;
    frame.set_atoms(make_water_atoms());
    auto uc = occ::crystal::triclinic_cell(a, b, c, alpha, beta, gamma);
    frame.set_unit_cell(uc);
    REQUIRE(writer.write_frame(frame));
    writer.finalise();
  }

  trajan::core::Frame frame;
  trajan::io::XYZHandler reader;
  reader.set_file_path(temp);
  reader.initialise(trajan::io::FileHandler::Mode::Read);
  REQUIRE(reader.read_frame(frame));
  reader.finalise();

  REQUIRE(frame.has_unit_cell());
  const auto &uc = frame.unit_cell().value();
  REQUIRE_THAT(uc.a(), Catch::Matchers::WithinAbs(a, 1e-4));
  REQUIRE_THAT(uc.b(), Catch::Matchers::WithinAbs(b, 1e-4));
  REQUIRE_THAT(uc.c(), Catch::Matchers::WithinAbs(c, 1e-4));
  REQUIRE_THAT(occ::units::degrees(uc.alpha()),
               Catch::Matchers::WithinAbs(80.0, 1e-3));
  REQUIRE_THAT(occ::units::degrees(uc.beta()),
               Catch::Matchers::WithinAbs(85.0, 1e-3));
  REQUIRE_THAT(occ::units::degrees(uc.gamma()),
               Catch::Matchers::WithinAbs(70.0, 1e-3));

  fs::remove(temp);
}

TEST_CASE("XYZ multi-frame round-trip", "[unit][io][xyz]") {
  fs::path temp = "/tmp/test_xyz_multiframe.xyz";

  {
    trajan::io::XYZHandler writer;
    writer.set_file_path(temp);
    writer.initialise(trajan::io::FileHandler::Mode::Write);

    for (int f = 0; f < 3; ++f) {
      trajan::core::Frame frame;
      auto atoms = make_water_atoms();
      for (auto &atom : atoms) {
        atom.x += static_cast<double>(f);
      }
      frame.set_atoms(atoms);
      REQUIRE(writer.write_frame(frame));
    }
    writer.finalise();
  }

  std::vector<std::vector<trajan::core::Atom>> all_atoms;
  {
    trajan::io::XYZHandler reader;
    reader.set_file_path(temp);
    reader.initialise(trajan::io::FileHandler::Mode::Read);

    trajan::core::Frame frame;
    while (reader.read_frame(frame)) {
      all_atoms.push_back(frame.atoms());
    }
    reader.finalise();
  }

  REQUIRE(all_atoms.size() == 3);
  for (int f = 0; f < 3; ++f) {
    REQUIRE(all_atoms[f].size() == 3);
    REQUIRE_THAT(all_atoms[f][0].x,
                 Catch::Matchers::WithinAbs(0.0 + f, 1e-5));
    REQUIRE_THAT(all_atoms[f][1].x,
                 Catch::Matchers::WithinAbs(0.960 + f, 1e-5));
  }

  fs::remove(temp);
}

TEST_CASE("XYZ via Trajectory - single frame", "[unit][io][xyz][trajectory]") {
  fs::path temp_write = "/tmp/test_traj_xyz_write.xyz";
  fs::path temp_read  = "/tmp/test_traj_xyz_read.xyz";

  {
    trajan::io::XYZHandler writer;
    writer.set_file_path(temp_write);
    writer.initialise(trajan::io::FileHandler::Mode::Write);

    trajan::core::Frame frame;
    frame.set_atoms(make_water_atoms());
    REQUIRE(writer.write_frame(frame));
    writer.finalise();
  }

  double ref_x0, ref_x1;
  {
    Trajectory traj;
    traj.load_files({temp_write});
    REQUIRE(traj.next_frame());

    const auto &atoms = traj.atoms();
    REQUIRE(atoms.size() == 3);
    REQUIRE(atoms[0].element.symbol() == "O");
    REQUIRE(atoms[1].element.symbol() == "H");
    REQUIRE(atoms[2].element.symbol() == "H");
    REQUIRE_THAT(atoms[0].x, Catch::Matchers::WithinAbs(0.000, 1e-5));
    REQUIRE_THAT(atoms[1].x, Catch::Matchers::WithinAbs(0.960, 1e-5));
    ref_x0 = atoms[0].x;
    ref_x1 = atoms[1].x;

    traj.set_output_file(temp_read);
    traj.write_frame();
    // traj destructor flushes and closes the output file
  }

  Trajectory traj2;
  traj2.load_files({temp_read});
  REQUIRE(traj2.next_frame());
  const auto &atoms2 = traj2.atoms();
  REQUIRE(atoms2.size() == 3);
  REQUIRE(atoms2[0].element.symbol() == "O");
  REQUIRE_THAT(atoms2[0].x, Catch::Matchers::WithinAbs(ref_x0, 1e-5));
  REQUIRE_THAT(atoms2[1].x, Catch::Matchers::WithinAbs(ref_x1, 1e-5));

  fs::remove(temp_write);
  fs::remove(temp_read);
}

TEST_CASE("XYZ via Trajectory - multi-frame with unit cell",
          "[unit][io][xyz][trajectory]") {
  fs::path temp = "/tmp/test_traj_xyz_uc.xyz";
  constexpr int n_frames = 5;
  constexpr double a = 20.0, b = 20.0, c = 20.0;

  {
    trajan::io::XYZHandler writer;
    writer.set_file_path(temp);
    writer.initialise(trajan::io::FileHandler::Mode::Write);

    auto uc = occ::crystal::triclinic_cell(a, b, c,
                                           occ::units::radians(90.0),
                                           occ::units::radians(90.0),
                                           occ::units::radians(90.0));
    for (int f = 0; f < n_frames; ++f) {
      trajan::core::Frame frame;
      auto atoms = make_water_atoms();
      for (auto &atom : atoms) {
        atom.x += 0.1 * f;
      }
      frame.set_atoms(atoms);
      frame.set_unit_cell(uc);
      REQUIRE(writer.write_frame(frame));
    }
    writer.finalise();
  }

  Trajectory traj;
  traj.load_files({temp});

  std::vector<double> ox_positions;
  while (traj.next_frame()) {
    ox_positions.push_back(traj.atoms()[0].x);
    REQUIRE(traj.unit_cell().has_value());
    REQUIRE_THAT(traj.unit_cell().value().a(),
                 Catch::Matchers::WithinAbs(a, 1e-4));
  }

  REQUIRE(static_cast<int>(ox_positions.size()) == n_frames);
  for (int f = 0; f < n_frames; ++f) {
    REQUIRE_THAT(ox_positions[f],
                 Catch::Matchers::WithinAbs(0.1 * f, 1e-5));
  }

  fs::remove(temp);
}

TEST_CASE("XYZ file - truncated atom list is rejected", "[unit][io][xyz]") {
  fs::path temp = "/tmp/test_xyz_corrupt.xyz";
  {
    std::ofstream f(temp);
    f << "5\n";
    f << "comment\n";
    f << "C 0.0 0.0 0.0\n";
    f << "H 1.1 0.0 0.0\n";
    f << "H -0.4 1.0 0.0\n";
    // only 3 of 5 claimed atoms
  }

  trajan::io::XYZHandler reader;
  reader.set_file_path(temp);
  reader.initialise(trajan::io::FileHandler::Mode::Read);

  trajan::core::Frame frame;
  REQUIRE_FALSE(reader.read_frame(frame));
  reader.finalise();

  fs::remove(temp);
}

TEST_CASE("XYZ file - diverse element symbols preserved", "[unit][io][xyz]") {
  fs::path temp = "/tmp/test_xyz_elements.xyz";
  {
    std::ofstream f(temp);
    f << "5\n";
    f << "Properties=species:S:1:pos:R:3\n";
    f << "C  0.0 0.0 0.0\n";
    f << "N  1.3 0.0 0.0\n";
    f << "O  2.5 0.0 0.0\n";
    f << "Fe 3.8 0.0 0.0\n";
    f << "Cl 5.1 0.0 0.0\n";
  }

  trajan::core::Frame frame;
  trajan::io::XYZHandler reader;
  reader.set_file_path(temp);
  reader.initialise(trajan::io::FileHandler::Mode::Read);
  REQUIRE(reader.read_frame(frame));
  reader.finalise();

  REQUIRE(frame.num_atoms() == 5);
  const auto &atoms = frame.atoms();
  REQUIRE(atoms[0].element.symbol() == "C");
  REQUIRE(atoms[1].element.symbol() == "N");
  REQUIRE(atoms[2].element.symbol() == "O");
  REQUIRE(atoms[3].element.symbol() == "Fe");
  REQUIRE(atoms[4].element.symbol() == "Cl");
  REQUIRE_THAT(atoms[3].x, Catch::Matchers::WithinAbs(3.8, 1e-5));

  fs::remove(temp);
}

// ── XTC tests (self-contained, no external test data needed) ─────────────────
//
// XTC stores coordinates in nm as float32 with lossy compression.
// With default precision=1000, the rounding error is ~0.5/1000 nm = 0.0005 nm
// = 0.005 Å.  A tolerance of 0.02 Å covers this comfortably.
// Box vectors are stored as raw float32 (no compression), so their tolerance
// is tighter (~1e-4 Å for typical box sizes).

static constexpr double XTC_COORD_TOL = 0.02;   // Angstroms
static constexpr double XTC_BOX_TOL   = 1e-3;   // Angstroms

// Helper: write a temp PDB so Trajectory has topology when loading XTC.
static fs::path write_topology_pdb(const fs::path &path,
                                   const std::vector<trajan::core::Atom> &atoms) {
  trajan::io::PDBHandler writer;
  writer.set_file_path(path);
  writer.initialise(trajan::io::FileHandler::Mode::Write);
  trajan::core::Frame frame;
  frame.set_atoms(atoms);
  writer.write_frame(frame);
  writer.finalise();
  return path;
}

TEST_CASE("XTC write/read round-trip - single frame, no unit cell",
          "[unit][io][xtc]") {
  fs::path temp = "/tmp/test_xtc_single.xtc";
  auto ref_atoms = make_water_atoms();

  // Write
  {
    trajan::io::XTCHandler writer;
    writer.set_file_path(temp);
    writer.initialise(trajan::io::FileHandler::Mode::Write);

    trajan::core::Frame frame;
    frame.set_atoms(ref_atoms);
    REQUIRE(writer.write_frame(frame));
    writer.finalise();
  }

  // Read back: pre-populate frame with atoms so update_atom_position works
  trajan::core::Frame frame;
  frame.set_atoms(ref_atoms);

  trajan::io::XTCHandler reader;
  reader.set_file_path(temp);
  reader.initialise(trajan::io::FileHandler::Mode::Read);
  REQUIRE(reader.read_frame(frame));

  // Positions should survive the nm round-trip within XTC tolerance
  const auto &atoms = frame.atoms();
  REQUIRE(atoms.size() == 3);
  REQUIRE_THAT(atoms[0].x, Catch::Matchers::WithinAbs(ref_atoms[0].x, XTC_COORD_TOL));
  REQUIRE_THAT(atoms[0].y, Catch::Matchers::WithinAbs(ref_atoms[0].y, XTC_COORD_TOL));
  REQUIRE_THAT(atoms[0].z, Catch::Matchers::WithinAbs(ref_atoms[0].z, XTC_COORD_TOL));
  REQUIRE_THAT(atoms[1].x, Catch::Matchers::WithinAbs(ref_atoms[1].x, XTC_COORD_TOL));
  REQUIRE_THAT(atoms[2].x, Catch::Matchers::WithinAbs(ref_atoms[2].x, XTC_COORD_TOL));
  REQUIRE_THAT(atoms[2].y, Catch::Matchers::WithinAbs(ref_atoms[2].y, XTC_COORD_TOL));

  // Second read returns false (only one frame written)
  trajan::core::Frame frame2;
  frame2.set_atoms(ref_atoms);
  REQUIRE_FALSE(reader.read_frame(frame2));
  reader.finalise();

  fs::remove(temp);
}

TEST_CASE("XTC write/read round-trip - single frame, cubic unit cell",
          "[unit][io][xtc]") {
  fs::path temp = "/tmp/test_xtc_uc.xtc";
  auto ref_atoms = make_water_atoms();
  constexpr double box_a = 20.0; // Angstroms

  {
    trajan::io::XTCHandler writer;
    writer.set_file_path(temp);
    writer.initialise(trajan::io::FileHandler::Mode::Write);

    trajan::core::Frame frame;
    frame.set_atoms(ref_atoms);
    auto uc = occ::crystal::triclinic_cell(box_a, box_a, box_a,
                                           occ::units::radians(90.0),
                                           occ::units::radians(90.0),
                                           occ::units::radians(90.0));
    frame.set_unit_cell(uc);
    REQUIRE(writer.write_frame(frame));
    writer.finalise();
  }

  trajan::core::Frame frame;
  frame.set_atoms(ref_atoms);

  trajan::io::XTCHandler reader;
  reader.set_file_path(temp);
  reader.initialise(trajan::io::FileHandler::Mode::Read);
  REQUIRE(reader.read_frame(frame));
  reader.finalise();

  REQUIRE(frame.has_unit_cell());
  const auto &uc = frame.unit_cell().value();
  REQUIRE_THAT(uc.a(), Catch::Matchers::WithinAbs(box_a, XTC_BOX_TOL));
  REQUIRE_THAT(uc.b(), Catch::Matchers::WithinAbs(box_a, XTC_BOX_TOL));
  REQUIRE_THAT(uc.c(), Catch::Matchers::WithinAbs(box_a, XTC_BOX_TOL));
  REQUIRE_THAT(occ::units::degrees(uc.alpha()),
               Catch::Matchers::WithinAbs(90.0, 1e-3));
  REQUIRE_THAT(occ::units::degrees(uc.beta()),
               Catch::Matchers::WithinAbs(90.0, 1e-3));
  REQUIRE_THAT(occ::units::degrees(uc.gamma()),
               Catch::Matchers::WithinAbs(90.0, 1e-3));

  fs::remove(temp);
}

TEST_CASE("XTC write/read round-trip - triclinic unit cell",
          "[unit][io][xtc]") {
  fs::path temp = "/tmp/test_xtc_triclinic.xtc";
  auto ref_atoms = make_water_atoms();
  constexpr double a = 15.0, b = 16.0, c = 17.0;
  const double alpha = occ::units::radians(80.0);
  const double beta  = occ::units::radians(85.0);
  const double gamma = occ::units::radians(75.0);

  {
    trajan::io::XTCHandler writer;
    writer.set_file_path(temp);
    writer.initialise(trajan::io::FileHandler::Mode::Write);

    trajan::core::Frame frame;
    frame.set_atoms(ref_atoms);
    auto uc = occ::crystal::triclinic_cell(a, b, c, alpha, beta, gamma);
    frame.set_unit_cell(uc);
    REQUIRE(writer.write_frame(frame));
    writer.finalise();
  }

  trajan::core::Frame frame;
  frame.set_atoms(ref_atoms);

  trajan::io::XTCHandler reader;
  reader.set_file_path(temp);
  reader.initialise(trajan::io::FileHandler::Mode::Read);
  REQUIRE(reader.read_frame(frame));
  reader.finalise();

  REQUIRE(frame.has_unit_cell());
  const auto &uc = frame.unit_cell().value();
  REQUIRE_THAT(uc.a(), Catch::Matchers::WithinAbs(a, XTC_BOX_TOL));
  REQUIRE_THAT(uc.b(), Catch::Matchers::WithinAbs(b, XTC_BOX_TOL));
  REQUIRE_THAT(uc.c(), Catch::Matchers::WithinAbs(c, XTC_BOX_TOL));
  REQUIRE_THAT(occ::units::degrees(uc.alpha()),
               Catch::Matchers::WithinAbs(80.0, 1e-2));
  REQUIRE_THAT(occ::units::degrees(uc.beta()),
               Catch::Matchers::WithinAbs(85.0, 1e-2));
  REQUIRE_THAT(occ::units::degrees(uc.gamma()),
               Catch::Matchers::WithinAbs(75.0, 1e-2));

  fs::remove(temp);
}

TEST_CASE("XTC multi-frame write/read round-trip", "[unit][io][xtc]") {
  fs::path temp = "/tmp/test_xtc_multi.xtc";
  constexpr int n_frames = 4;
  auto ref_atoms = make_water_atoms();

  // Write n_frames frames with linearly displaced O atom
  {
    trajan::io::XTCHandler writer;
    writer.set_file_path(temp);
    writer.initialise(trajan::io::FileHandler::Mode::Write);

    for (int f = 0; f < n_frames; ++f) {
      trajan::core::Frame frame;
      auto atoms = make_water_atoms();
      for (auto &atom : atoms) {
        atom.x += 0.5 * f; // 0.5 Å shift per frame
        atom.y += 0.3 * f;
      }
      frame.set_atoms(atoms);
      REQUIRE(writer.write_frame(frame));
    }
    writer.finalise();
  }

  // Read all frames back
  std::vector<double> ox_x, ox_y;
  {
    trajan::io::XTCHandler reader;
    reader.set_file_path(temp);
    reader.initialise(trajan::io::FileHandler::Mode::Read);

    trajan::core::Frame frame;
    frame.set_atoms(ref_atoms); // pre-load atom metadata
    while (reader.read_frame(frame)) {
      ox_x.push_back(frame.atoms()[0].x);
      ox_y.push_back(frame.atoms()[0].y);
    }
    reader.finalise();
  }

  REQUIRE(static_cast<int>(ox_x.size()) == n_frames);
  for (int f = 0; f < n_frames; ++f) {
    REQUIRE_THAT(ox_x[f],
                 Catch::Matchers::WithinAbs(0.0 + 0.5 * f, XTC_COORD_TOL));
    REQUIRE_THAT(ox_y[f],
                 Catch::Matchers::WithinAbs(0.0 + 0.3 * f, XTC_COORD_TOL));
  }

  fs::remove(temp);
}

TEST_CASE("XTC via Trajectory - PDB topology + XTC coordinates",
          "[unit][io][xtc][trajectory]") {
  fs::path temp_pdb = "/tmp/test_xtc_traj.pdb";
  fs::path temp_xtc = "/tmp/test_xtc_traj.xtc";
  constexpr int n_frames = 3;
  auto ref_atoms = make_water_atoms();

  write_topology_pdb(temp_pdb, ref_atoms);

  // Write XTC frames with known positions
  {
    trajan::io::XTCHandler writer;
    writer.set_file_path(temp_xtc);
    writer.initialise(trajan::io::FileHandler::Mode::Write);

    auto uc = occ::crystal::triclinic_cell(25.0, 25.0, 25.0,
                                           occ::units::radians(90.0),
                                           occ::units::radians(90.0),
                                           occ::units::radians(90.0));
    for (int f = 0; f < n_frames; ++f) {
      trajan::core::Frame frame;
      auto atoms = make_water_atoms();
      for (auto &atom : atoms) {
        atom.x += 1.0 * f;
      }
      frame.set_atoms(atoms);
      frame.set_unit_cell(uc);
      REQUIRE(writer.write_frame(frame));
    }
    writer.finalise();
  }

  // Load via Trajectory (PDB provides topology, XTC provides coordinates)
  Trajectory traj;
  traj.load_files({temp_pdb, temp_xtc});

  // First next_frame() reads PDB topology (atom types/connectivity, no unit cell)
  REQUIRE(traj.next_frame());
  REQUIRE(traj.atoms().size() == 3);

  // Subsequent next_frame() calls read XTC coordinate frames
  std::vector<double> ox_x;
  while (traj.next_frame()) {
    ox_x.push_back(traj.atoms()[0].x);
    REQUIRE(traj.atoms().size() == 3);
    REQUIRE(traj.unit_cell().has_value());
    REQUIRE_THAT(traj.unit_cell().value().a(),
                 Catch::Matchers::WithinAbs(25.0, XTC_BOX_TOL));
  }

  REQUIRE(static_cast<int>(ox_x.size()) == n_frames);
  for (int f = 0; f < n_frames; ++f) {
    REQUIRE_THAT(ox_x[f],
                 Catch::Matchers::WithinAbs(0.0 + 1.0 * f, XTC_COORD_TOL));
  }

  fs::remove(temp_pdb);
  fs::remove(temp_xtc);
}

TEST_CASE("XTC via Trajectory - write output from PDB+XTC, re-read",
          "[unit][io][xtc][trajectory]") {
  fs::path temp_pdb  = "/tmp/test_xtc_rewrite.pdb";
  fs::path temp_xtc  = "/tmp/test_xtc_rewrite.xtc";
  fs::path temp_xtc2 = "/tmp/test_xtc_rewrite2.xtc";
  auto ref_atoms = make_water_atoms();

  write_topology_pdb(temp_pdb, ref_atoms);

  // Write a 2-frame XTC
  {
    trajan::io::XTCHandler writer;
    writer.set_file_path(temp_xtc);
    writer.initialise(trajan::io::FileHandler::Mode::Write);
    for (int f = 0; f < 2; ++f) {
      trajan::core::Frame frame;
      auto atoms = make_water_atoms();
      for (auto &atom : atoms) atom.x += 2.0 * f;
      frame.set_atoms(atoms);
      REQUIRE(writer.write_frame(frame));
    }
    writer.finalise();
  }

  // Read via Trajectory, write each XTC frame out as a second XTC
  // (first next_frame() loads PDB topology and is not written)
  {
    Trajectory traj;
    traj.load_files({temp_pdb, temp_xtc});
    traj.set_output_file(temp_xtc2);
    REQUIRE(traj.next_frame()); // load PDB topology, skip writing
    while (traj.next_frame()) {
      traj.write_frame();
    }
    // traj destructor flushes and closes temp_xtc2
  }

  // Re-read the rewritten XTC
  std::vector<double> ox_x;
  {
    trajan::io::XTCHandler reader;
    reader.set_file_path(temp_xtc2);
    reader.initialise(trajan::io::FileHandler::Mode::Read);

    trajan::core::Frame frame;
    frame.set_atoms(ref_atoms);
    while (reader.read_frame(frame)) {
      ox_x.push_back(frame.atoms()[0].x);
    }
    reader.finalise();
  }

  REQUIRE(ox_x.size() == 2);
  REQUIRE_THAT(ox_x[0], Catch::Matchers::WithinAbs(0.0, XTC_COORD_TOL));
  REQUIRE_THAT(ox_x[1], Catch::Matchers::WithinAbs(2.0, XTC_COORD_TOL));

  fs::remove(temp_pdb);
  fs::remove(temp_xtc);
  fs::remove(temp_xtc2);
}

TEST_CASE_METHOD(TestFixture, "XTC Read - external test file",
                 "[file][io][xtc]") {
  // Uses the md.xtc from the libxtc examples if available via --data-path,
  // otherwise skips.  Checks basic frame iteration and atom count consistency.
  fs::path xtc_file = this->get_test_file("md.xtc");
  fs::path pdb_file = this->get_test_file("md.pdb");

  // Need both topology and trajectory
  if (!fs::exists(pdb_file)) {
    SKIP("md.pdb not found alongside md.xtc — skipping external XTC test");
  }

  Trajectory traj;
  traj.load_files({pdb_file, xtc_file});

  size_t frame_count = 0;
  size_t natoms = 0;
  while (traj.next_frame()) {
    if (frame_count == 0) {
      natoms = traj.atoms().size();
      REQUIRE(natoms > 0);
    } else {
      REQUIRE(traj.atoms().size() == natoms);
    }
    ++frame_count;
  }
  REQUIRE(frame_count > 0);
}
