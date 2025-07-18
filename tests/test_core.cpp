#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <chrono>
#include <omp.h>
#include <random>
#include <trajan/core/element.h>
#include <trajan/core/linear_algebra.h>
#include <trajan/core/molecule.h>
#include <trajan/core/neigh.h>
#include <trajan/core/topology.h>
#include <trajan/core/trajectory.h>
#include <trajan/core/unit_cell.h>
#include <trajan/core/units.h>
#include <trajan/io/pdb.h>

namespace fs = std::filesystem;
#include <trajan/io/selection.h>
#include <vector>

using trajan::Mat3N;
using trajan::Vec3;
using trajan::core::Atom;
using trajan::core::Bond;
using trajan::core::BondGraph;
using trajan::core::Cell;
using trajan::core::CellList;
using trajan::core::DihedralType;
using trajan::core::Element;
using trajan::core::Entity;
using trajan::core::Frame;
using trajan::core::Molecule;
using trajan::core::Topology;
using trajan::core::Trajectory;
using trajan::core::UnitCell;
using trajan::io::SelectionCriteria;
using trajan::io::SelectionParser;
using Entities = std::vector<trajan::core::EntityType>;

namespace units = trajan::units;

fs::path CURRENT = __FILE__;
fs::path EXAMPLES_DIR = CURRENT.parent_path().parent_path() / "examples";

TEST_CASE("Element constructor with exact matching", "[element]") {
  SECTION("Exact matches are found correctly") {
    auto [input, expected_symbol, expected_number] =
        GENERATE(table<std::string, std::string, int>({{"H", "H", 1},
                                                       {"He", "He", 2},
                                                       {"Fe", "Fe", 26},
                                                       {"Au", "Au", 79},
                                                       {"U", "U", 92}}));

    Element e(input, true);
    REQUIRE(e.symbol() == expected_symbol);
    REQUIRE(e.atomic_number() == expected_number);
  }

  SECTION("Case sensitivity for exact matches") {
    auto [input, expected_symbol] = GENERATE(table<std::string, std::string>(
        {{"h", "Xx"}, {"HE", "Xx"}, {"fe", "Xx"}, {"AU", "Xx"}}));

    Element e(input, true);
    REQUIRE(e.symbol() == expected_symbol);
  }
}

TEST_CASE("Element constructor with partial matching", "[element]") {
  SECTION("Partial matches find longest valid symbol") {
    auto [input, expected_symbol] =
        GENERATE(table<std::string, std::string>({{"Hg", "Hg"},
                                                  {"Hge", "Hg"},
                                                  {"Her", "He"},
                                                  {"Fer", "Fe"},
                                                  {"Feat", "Fe"}}));

    Element e(input, false);
    REQUIRE(e.symbol() == expected_symbol);
  }

  SECTION("No match returns dummy element") {
    auto input = GENERATE("Q", "J", "X", "Qq", "Jk", "Xyz");

    Element e(input, false);
    REQUIRE(e.symbol() == "Xx");
    REQUIRE(e.atomic_number() == 0);
  }
}

TEST_CASE("Element constructor with edge cases", "[element]") {
  SECTION("Empty string returns dummy element") {
    Element e("", false);
    REQUIRE(e.symbol() == "Xx");
  }

  SECTION("Single character matches") {
    auto [input, expected_symbol] = GENERATE(table<std::string, std::string>(
        {{"H", "H"}, {"N", "N"}, {"O", "O"}, {"F", "F"}, {"U", "U"}}));

    Element e(input, false);
    REQUIRE(e.symbol() == expected_symbol);
  }
}

TEST_CASE("CellList Basic Construction", "[cell_list]") {
  auto unit_cell = trajan::core::cubic_cell(10.0);
  double cutoff = 2.0;

  SECTION("Construction with single thread") {
    REQUIRE_NOTHROW(CellList(unit_cell, cutoff, 1));
  }

  SECTION("Construction with multiple threads") {
    REQUIRE_NOTHROW(CellList(unit_cell, cutoff, 4));
  }

  SECTION("Construction with invalid thread count") {
    REQUIRE_NOTHROW(CellList(unit_cell, cutoff, -1));
  }
}

TEST_CASE("Cell Adding and Retrieving Atoms", "[cell]") {
  Cell cell;
  Vec3 pos(1.0, 1.0, 1.0);
  Atom atom(pos, 0, 0);

  SECTION("Adding single atom") {
    cell.add_entity(atom);
    REQUIRE(cell.get_entities().size() == 1);
    REQUIRE(cell.get_entities()[0].position() == pos);
  }

  SECTION("Adding multiple atoms") {
    const int num_atoms = 10;
    std::vector<Atom> atoms;
    for (int i = 0; i < num_atoms; ++i) {
      atoms.emplace_back(Vec3(i, i, i), 0, i);
      cell.add_entity(atoms.back());
    }
    REQUIRE(cell.get_entities().size() == num_atoms);
  }
}

// helper function to measure execution time
template <typename Func> double measure_execution_time(Func &&func) {
  auto start = std::chrono::high_resolution_clock::now();
  func();
  auto end = std::chrono::high_resolution_clock::now();
  return std::chrono::duration<double>(end - start).count();
}

TEST_CASE("CellList vs Double Loop Comparison", "[cell_list]") {
  const double box_size = 50.0;
  UnitCell unit_cell = trajan::core::cubic_cell(box_size);
  const int num_atoms = 3000;
  const double cutoff = 5.0;
  const int num_threads = 1;
  std::vector<Atom> atoms;
  atoms.reserve(num_atoms);
  std::random_device rd;
  std::mt19937 gen(42);

  SECTION("Comparing pair counting between methods inside box only atoms") {
    std::uniform_real_distribution<> dis(0, box_size);

    // create atoms with random positions
    for (int i = 0; i < num_atoms; ++i) {
      Vec3 position(dis(gen), dis(gen), dis(gen));
      Atom atom(position, 0, i);
      atoms.push_back(atom);
    }

    std::atomic<size_t> cell_list_pairs{0};
    std::atomic<size_t> verlet_list_pairs{0};

    trajan::core::NeighbourList neighbour_list(unit_cell, cutoff, num_threads);
    double cell_list_time = measure_execution_time([&]() {
      neighbour_list.update(atoms);
      neighbour_list.iterate_neighbours([&](const Entity &e1, const Entity &e2,
                                            double rsq) { cell_list_pairs++; });
    });

    neighbour_list =
        trajan::core::NeighbourList(unit_cell, cutoff, num_threads,
                                    trajan::core::NeighbourList::Type::Verlet);
    double verlet_list_time = measure_execution_time([&]() {
      neighbour_list.update(atoms);
      neighbour_list.iterate_neighbours(
          [&](const Entity &e1, const Entity &e2, double rsq) {
            verlet_list_pairs++;
          });
    });

    INFO("System size: " << num_atoms << " atoms");
    INFO("Box size: " << box_size << "x" << box_size << "x" << box_size);
    INFO("Cutoff distance: " << cutoff);
    INFO("Number of threads: " << num_threads);
    double volume = box_size * box_size * box_size;
    double expected_pairs_per_atom =
        (4.0 / 3.0) * units::PI * std::pow(cutoff, 3) / volume * num_atoms;
    double total_expected_pairs = expected_pairs_per_atom * num_atoms / 2;
    INFO("Expected approximate number of pairs: ~" << total_expected_pairs);
    INFO("Pairs found (cell list): " << cell_list_pairs);
    INFO("Pairs found (double loop): " << verlet_list_pairs);
    REQUIRE(cell_list_pairs == verlet_list_pairs);
    REQUIRE(static_cast<double>(cell_list_pairs) > total_expected_pairs * 0.8);
    REQUIRE(static_cast<double>(cell_list_pairs) < total_expected_pairs * 1.2);
    INFO("Cell list execution time: " << cell_list_time << " seconds");
    INFO("Verlet list execution time: " << verlet_list_time << " seconds");
    INFO("Speed-up factor: " << verlet_list_time / cell_list_time);
    REQUIRE(verlet_list_time > cell_list_time);
  }
  SECTION("Comparing pair counting between methods outside box only atoms") {
    std::uniform_real_distribution<> dis(-5, box_size - 5);

    atoms.clear();
    atoms.reserve(num_atoms);
    // create atoms with random positions
    for (int i = 0; i < num_atoms; ++i) {
      Vec3 position(dis(gen), dis(gen), dis(gen));
      Atom atom(position, 0, i);
      atoms.push_back(atom);
    }
    std::atomic<size_t> cell_list_pairs{0};
    std::atomic<size_t> verlet_list_pairs{0};

    trajan::core::NeighbourList neighbour_list(unit_cell, cutoff, num_threads);
    double cell_list_time = measure_execution_time([&]() {
      neighbour_list.update(atoms);
      neighbour_list.iterate_neighbours([&](const Entity &e1, const Entity &e2,
                                            double rsq) { cell_list_pairs++; });
    });

    neighbour_list =
        trajan::core::NeighbourList(unit_cell, cutoff, num_threads,
                                    trajan::core::NeighbourList::Type::Verlet);
    double verlet_list_time = measure_execution_time([&]() {
      neighbour_list.update(atoms);
      neighbour_list.iterate_neighbours(
          [&](const Entity &e1, const Entity &e2, double rsq) {
            verlet_list_pairs++;
          });
    });

    INFO("System size: " << num_atoms << " atoms");
    INFO("Box size: " << box_size << "x" << box_size << "x" << box_size);
    INFO("Cutoff distance: " << cutoff);
    INFO("Number of threads: " << num_threads);
    double volume = box_size * box_size * box_size;
    double expected_pairs_per_atom =
        (4.0 / 3.0) * units::PI * std::pow(cutoff, 3) / volume * num_atoms;
    double total_expected_pairs = expected_pairs_per_atom * num_atoms / 2;
    INFO("Expected approximate number of pairs: ~" << total_expected_pairs);
    INFO("Pairs found (cell list): " << cell_list_pairs);
    INFO("Pairs found (double loop): " << verlet_list_pairs);
    REQUIRE(cell_list_pairs == verlet_list_pairs);
    REQUIRE(static_cast<double>(cell_list_pairs) > total_expected_pairs * 0.8);
    REQUIRE(static_cast<double>(cell_list_pairs) < total_expected_pairs * 1.2);
    INFO("Cell list execution time: " << cell_list_time << " seconds");
    INFO("Verlet list execution time: " << verlet_list_time << " seconds");
    INFO("Speed-up factor: " << verlet_list_time / cell_list_time);
    REQUIRE(verlet_list_time > cell_list_time);
  }
}
TEST_CASE("CellList vs Double Loop Comparison using entities", "[cell_list]") {
  const double box_size = 50.0;
  const int num_atoms = 3000;
  const double cutoff = 5.0;
  const int num_threads = 1;
  std::random_device rd;
  std::mt19937 gen(42);
  std::uniform_real_distribution<> dis(0, box_size);

  UnitCell unit_cell = trajan::core::cubic_cell(box_size);

  // create atoms with random positions
  Entities entities;
  entities.reserve(num_atoms);
  for (int i = 0; i < num_atoms; ++i) {
    Vec3 position(dis(gen), dis(gen), dis(gen));
    Atom atom(position, 0, i);
    entities.push_back(atom);
  }

  SECTION("Comparing pair counting between methods using entities") {
    std::atomic<size_t> cell_list_pairs{0};
    std::atomic<size_t> verlet_list_pairs{0};

    trajan::core::NeighbourList neighbour_list(unit_cell, cutoff, num_threads);
    double cell_list_time = measure_execution_time([&]() {
      neighbour_list.update(entities);
      neighbour_list.iterate_neighbours([&](const Entity &e1, const Entity &e2,
                                            double rsq) { cell_list_pairs++; });
    });

    neighbour_list =
        trajan::core::NeighbourList(unit_cell, cutoff, num_threads,
                                    trajan::core::NeighbourList::Type::Verlet);
    double verlet_list_time = measure_execution_time([&]() {
      neighbour_list.update(entities);
      neighbour_list.iterate_neighbours(
          [&](const Entity &e1, const Entity &e2, double rsq) {
            verlet_list_pairs++;
          });
    });

    INFO("System size: " << num_atoms << " atoms");
    INFO("Box size: " << box_size << "x" << box_size << "x" << box_size);
    INFO("Cutoff distance: " << cutoff);
    INFO("Number of threads: " << num_threads);
    double volume = box_size * box_size * box_size;
    double expected_pairs_per_atom =
        (4.0 / 3.0) * units::PI * std::pow(cutoff, 3) / volume * num_atoms;
    double total_expected_pairs = expected_pairs_per_atom * num_atoms / 2;
    INFO("Expected approximate number of pairs: ~" << total_expected_pairs);
    INFO("Pairs found (cell list): " << cell_list_pairs);
    INFO("Pairs found (double loop): " << verlet_list_pairs);
    REQUIRE(cell_list_pairs == verlet_list_pairs);
    REQUIRE(static_cast<double>(cell_list_pairs) > total_expected_pairs * 0.8);
    REQUIRE(static_cast<double>(cell_list_pairs) < total_expected_pairs * 1.2);
    INFO("Cell list execution time: " << cell_list_time << " seconds");
    INFO("Verlet list execution time: " << verlet_list_time << " seconds");
    INFO("Speed-up factor: " << verlet_list_time / cell_list_time);
    REQUIRE(verlet_list_time > cell_list_time);
  }
}
TEST_CASE("CellList vs Double Loop Comparison using entities in Trajectory",
          "[cell_list]") {
  const double box_size = 50.0;
  const int num_atoms = 3000;
  const double cutoff = 5.0;
  const int num_threads = 1;
  std::random_device rd;
  std::mt19937 gen(42);
  std::uniform_real_distribution<> dis(0, box_size);

  Trajectory trajectory;
  Frame &frame = trajectory.frame();
  UnitCell unit_cell = trajan::core::cubic_cell(box_size);
  frame.set_uc(unit_cell);

  // create atoms with random positions
  std::vector<Atom> atoms;
  atoms.reserve(num_atoms);
  for (int i = 0; i < num_atoms; ++i) {
    Vec3 position(dis(gen), dis(gen), dis(gen));
    Atom atom(position, 0, i);
    atoms.push_back(atom);
  }
  frame.set_atoms(atoms);

  SECTION("Comparing pair counting between methods using entities") {
    std::atomic<size_t> cell_list_pairs{0};
    std::atomic<size_t> verlet_list_pairs{0};

    trajan::core::NeighbourList neighbour_list(unit_cell, cutoff, num_threads);
    double cell_list_time = measure_execution_time([&]() {
      neighbour_list.update(trajectory.atoms());
      neighbour_list.iterate_neighbours([&](const Entity &e1, const Entity &e2,
                                            double rsq) { cell_list_pairs++; });
    });

    neighbour_list =
        trajan::core::NeighbourList(unit_cell, cutoff, num_threads,
                                    trajan::core::NeighbourList::Type::Verlet);
    double verlet_list_time = measure_execution_time([&]() {
      neighbour_list.update(trajectory.atoms());
      neighbour_list.iterate_neighbours(
          [&](const Entity &e1, const Entity &e2, double rsq) {
            verlet_list_pairs++;
          });
    });

    INFO("System size: " << num_atoms << " atoms");
    INFO("Box size: " << box_size << "x" << box_size << "x" << box_size);
    INFO("Cutoff distance: " << cutoff);
    INFO("Number of threads: " << num_threads);
    double volume = box_size * box_size * box_size;
    double expected_pairs_per_atom =
        (4.0 / 3.0) * units::PI * std::pow(cutoff, 3) / volume * num_atoms;
    double total_expected_pairs = expected_pairs_per_atom * num_atoms / 2;
    INFO("Expected approximate number of pairs: ~" << total_expected_pairs);
    INFO("Pairs found (cell list): " << cell_list_pairs);
    INFO("Pairs found (double loop): " << verlet_list_pairs);
    REQUIRE(cell_list_pairs == verlet_list_pairs);
    REQUIRE(static_cast<double>(cell_list_pairs) > total_expected_pairs * 0.8);
    REQUIRE(static_cast<double>(cell_list_pairs) < total_expected_pairs * 1.2);
    INFO("Cell list execution time: " << cell_list_time << " seconds");
    INFO("Verlet list execution time: " << verlet_list_time << " seconds");
    INFO("Speed-up factor: " << verlet_list_time / cell_list_time);
    REQUIRE(verlet_list_time > cell_list_time);
  }
}

TEST_CASE("CellList vs Double Loop Comparison inside Trajectory with selection",
          "[cell_list]") {
  const double box_size = 50.0;
  const int num_atoms = 3000;
  const double cutoff = 5.0;
  const int num_threads = 1;
  std::random_device rd;
  std::mt19937 gen(42);
  std::uniform_real_distribution<> dis(0, box_size);

  Trajectory trajectory;
  Frame &frame = trajectory.frame();
  UnitCell unit_cell = trajan::core::cubic_cell(box_size);
  frame.set_uc(unit_cell);
  std::vector<Atom> atoms;
  for (int i = 0; i < num_atoms; ++i) {
    Vec3 position(dis(gen), dis(gen), dis(gen));
    Atom atom;
    atom.x = position.x();
    atom.y = position.y();
    atom.z = position.z();
    atom.index = i;
    if (i % 3 == 0) {
      atom.type = "O2";
    } else {
      atom.type = "H2";
    }
    atoms.push_back(atom);
  }
  frame.set_atoms(atoms);
  std::string sel = "aO2";
  std::optional<SelectionCriteria> result = SelectionParser::parse(sel);
  REQUIRE(result);
  SelectionCriteria sc = *result;
  Entities entities = trajectory.get_entities(sc);
  size_t entities_size = entities.size();

  SECTION("Comparing pair counting between methods") {
    std::atomic<size_t> cell_list_pairs{0};
    std::atomic<size_t> verlet_list_pairs{0};

    trajan::core::NeighbourList neighbour_list(unit_cell, cutoff, num_threads);
    double cell_list_time = measure_execution_time([&]() {
      neighbour_list.update(entities);
      neighbour_list.iterate_neighbours([&](const Entity &e1, const Entity &e2,
                                            double rsq) { cell_list_pairs++; });
    });

    neighbour_list =
        trajan::core::NeighbourList(unit_cell, cutoff, num_threads,
                                    trajan::core::NeighbourList::Type::Verlet);
    double verlet_list_time = measure_execution_time([&]() {
      neighbour_list.update(entities);
      neighbour_list.iterate_neighbours(
          [&](const Entity &e1, const Entity &e2, double rsq) {
            verlet_list_pairs++;
          });
    });

    INFO("System size: " << num_atoms << " atoms");
    INFO("Number of entities: " << entities_size);
    INFO("Box size: " << box_size << "x" << box_size << "x" << box_size);
    INFO("Cutoff distance: " << cutoff);
    INFO("Number of threads: " << num_threads);
    double volume = box_size * box_size * box_size;
    double expected_pairs_per_atom =
        (4.0 / 3.0) * units::PI * std::pow(cutoff, 3) / volume * entities_size;
    double total_expected_pairs = expected_pairs_per_atom * entities_size / 2;
    INFO("Expected approximate number of pairs: ~" << total_expected_pairs);
    INFO("Pairs found (cell list): " << cell_list_pairs);
    INFO("Pairs found (double loop): " << verlet_list_pairs);
    REQUIRE(cell_list_pairs == verlet_list_pairs);
    REQUIRE(static_cast<double>(cell_list_pairs) > total_expected_pairs * 0.8);
    REQUIRE(static_cast<double>(cell_list_pairs) < total_expected_pairs * 1.2);
    INFO("Cell list execution time: " << cell_list_time << " seconds");
    INFO("Verlet list execution time: " << verlet_list_time << " seconds");
    INFO("Speed-up factor: " << verlet_list_time / cell_list_time);
    REQUIRE(verlet_list_time > cell_list_time);
  }
}
TEST_CASE(
    "CellList vs Double Loop Comparison inside Trajectory with multi selection",
    "[cell_list]") {
  const double box_size = 25.0;
  const int num_atoms = 1000;
  const double cutoff = 6.0;
  const int num_threads = 1;
  Trajectory trajectory;
  Frame &frame = trajectory.frame();
  UnitCell unit_cell = trajan::core::cubic_cell(box_size);
  frame.set_uc(unit_cell);
  std::vector<Atom> atoms;

  int points_per_dim = static_cast<int>(std::ceil(std::cbrt(num_atoms)));

  double spacing = box_size / points_per_dim;

  for (int i = 0; i < num_atoms; ++i) {
    int x_idx = i % points_per_dim;
    int y_idx = (i / points_per_dim) % points_per_dim;
    int z_idx = i / (points_per_dim * points_per_dim);

    double x = (x_idx + 0.5) * spacing;
    double y = (y_idx + 0.5) * spacing;
    double z = (z_idx + 0.5) * spacing;

    Atom atom;
    atom.x = x;
    atom.y = y;
    atom.z = z;
    atom.index = i;
    atom.type = "O2";

    atoms.push_back(atom);
  }

  frame.set_atoms(atoms);
  std::string sel = "aO2";
  std::optional<SelectionCriteria> result1 = SelectionParser::parse(sel);
  REQUIRE(result1);
  SelectionCriteria sc1 = *result1;
  std::optional<SelectionCriteria> result2 = SelectionParser::parse(sel);
  REQUIRE(result2);
  SelectionCriteria sc2 = *result2;
  std::vector<Entities> all_entities = {trajectory.get_entities(sc1),
                                        trajectory.get_entities(sc2)};
  size_t entities_size = all_entities[0].size();

  SECTION("Comparing pair counting between methods") {
    std::atomic<size_t> cell_list_pairs{0};
    std::atomic<size_t> verlet_list_pairs{0};

    trajan::core::NeighbourList neighbour_list(unit_cell, cutoff, num_threads);
    double cell_list_time = measure_execution_time([&]() {
      neighbour_list.update(all_entities);
      neighbour_list.iterate_neighbours([&](const trajan::core::Entity &e1,
                                            const trajan::core::Entity &e2,
                                            double rsq) { cell_list_pairs++; });
    });

    neighbour_list =
        trajan::core::NeighbourList(unit_cell, cutoff, num_threads,
                                    trajan::core::NeighbourList::Type::Verlet);
    double verlet_list_time = measure_execution_time([&]() {
      neighbour_list.update(all_entities);
      neighbour_list.iterate_neighbours(
          [&](const trajan::core::Entity &e1, const trajan::core::Entity &e2,
              double rsq) { verlet_list_pairs++; });
    });

    INFO("System size: " << num_atoms << " atoms");
    INFO("Number of entities: " << entities_size);
    INFO("Box size: " << box_size << "x" << box_size << "x" << box_size);
    INFO("Cutoff distance: " << cutoff);
    INFO("Number of threads: " << num_threads);
    double volume = box_size * box_size * box_size;
    double expected_pairs_per_atom =
        (4.0 / 3.0) * units::PI * std::pow(cutoff, 3) / volume * entities_size;
    double total_expected_pairs = expected_pairs_per_atom * entities_size;
    INFO("Expected approximate number of pairs: ~" << total_expected_pairs);
    INFO("Pairs found (cell list): " << cell_list_pairs);
    INFO("Pairs found (double loop): " << verlet_list_pairs);
    REQUIRE(cell_list_pairs == verlet_list_pairs);
    REQUIRE(static_cast<double>(cell_list_pairs) > total_expected_pairs * 0.8);
    REQUIRE(static_cast<double>(cell_list_pairs) < total_expected_pairs * 1.2);
    INFO("Cell list execution time: " << cell_list_time << " seconds");
    INFO("Verlet list execution time: " << verlet_list_time << " seconds");
    INFO("Speed-up factor: " << verlet_list_time / cell_list_time);
    REQUIRE(verlet_list_time > cell_list_time);
  }
}

TEST_CASE("Molecule construction from atoms", "[molecule]") {
  SECTION("Empty molecule") {
    std::vector<Atom> atoms;
    Molecule molecule(atoms);

    REQUIRE(molecule.size() == 0);
  }

  SECTION("Single atom molecule") {
    std::vector<Atom> atoms = {
        Atom(Vec3{0.0, 0.0, 0.0}, Element("H", true), 0)};

    Molecule molecule(atoms);

    REQUIRE(molecule.size() == 1);
  }

  SECTION("Water molecule") {
    std::vector<Atom> atoms = {
        Atom(Vec3{0.0, 0.0, 0.0}, Element("O", true), 0),
        Atom(Vec3{0.8, 0.0, 0.0}, Element("H", true), 1),
        Atom(Vec3{-0.2, 0.8, 0.0}, Element("H", true), 2)};

    Molecule molecule(atoms);

    REQUIRE(molecule.size() == 3);
  }
}

TEST_CASE("MoleculeGraph bond detection", "[molecule][graph]") {
  SECTION("Bonded atoms") {
    // Typical C-O bond length
    Atom a1(Vec3{0.0, 0.0, 0.0}, Element("C", true), 0);
    Atom a2(Vec3{1.2, 0.0, 0.0}, Element("O", true), 1);

    auto bond_opt = a1.is_bonded(a2);

    REQUIRE(bond_opt.has_value());
    REQUIRE(bond_opt->bond_length == Catch::Approx(1.2));
  }

  SECTION("Non-bonded atoms") {
    // Too far for bonding
    Atom a1(Vec3{0.0, 0.0, 0.0}, Element("C", true), 0);
    Atom a2(Vec3{5.0, 0.0, 0.0}, Element("O", true), 1);

    auto bond_opt = a1.is_bonded(a2);

    REQUIRE_FALSE(bond_opt.has_value());
  }

  SECTION("Boundary case with tolerance") {
    // Slightly longer than typical C-O bond
    Atom a1(Vec3{0.0, 0.0, 0.0}, Element("C", true), 0);
    Atom a2(Vec3{1.8, 0.0, 0.0}, Element("O", true), 1);

    // With default tolerance
    auto bond_default = a1.is_bonded(a2);
    // With increased tolerance
    auto bond_increased = a1.is_bonded(a2, 0.6);

    REQUIRE_FALSE(bond_default.has_value());
    REQUIRE(bond_increased.has_value());
  }
}

std::vector<Atom> create_test_atoms() {
  std::vector<Atom> atoms;

  // Create a simple methane-like molecule (C-H4)
  atoms.emplace_back(Vec3{0.0, 0.0, 0.0}, "C",
                     0);                            // Carbon center
  atoms.emplace_back(Vec3{1.1, 0.0, 0.0}, "H", 1);  // H1
  atoms.emplace_back(Vec3{-1.1, 0.0, 0.0}, "H", 2); // H2
  atoms.emplace_back(Vec3{0.0, 1.1, 0.0}, "H", 3);  // H3
  atoms.emplace_back(Vec3{0.0, -1.1, 0.0}, "H", 4); // H4

  return atoms;
}

std::vector<Atom> create_linear_molecule() {
  std::vector<Atom> atoms;

  // Create a linear molecule H-C-C-H
  atoms.emplace_back(Vec3{0.0, 0.0, 0.0}, "H", 0);
  atoms.emplace_back(Vec3{1.1, 0.0, 0.0}, "C", 1);
  atoms.emplace_back(Vec3{2.6, 0.0, 0.0}, "C", 2);
  atoms.emplace_back(Vec3{3.7, 0.0, 0.0}, "H", 3);

  return atoms;
}

std::vector<Atom> create_disconnected_atoms() {
  std::vector<Atom> atoms;

  // Create two separate H2 molecules
  atoms.emplace_back(Vec3{0.0, 0.0, 0.0}, "H", 0);
  atoms.emplace_back(Vec3{0.8, 0.0, 0.0}, "H", 1);
  atoms.emplace_back(Vec3{10.0, 0.0, 0.0}, "H", 2);
  atoms.emplace_back(Vec3{10.8, 0.0, 0.0}, "H", 3);

  return atoms;
}

TEST_CASE("Topology Construction", "[topology][construction]") {
  SECTION("Default constructor") {
    Topology topology;
    REQUIRE(topology.num_bonds() == 0);
    REQUIRE(topology.num_angles() == 0);
    REQUIRE(topology.num_dihedrals() == 0);
  }

  SECTION("Construction from atoms with automatic bond detection") {
    auto atoms = create_test_atoms();
    Topology topology(atoms);

    // Should detect C-H bonds automatically
    REQUIRE(topology.num_bonds() > 0);
    REQUIRE(topology.num_angles() > 0);
  }

  SECTION("Construction from BondGraph") {
    auto atoms = create_test_atoms();
    BondGraph bond_graph(atoms);

    // Manually add some bonds
    Bond bond(1.1);
    bond_graph.add_edge(0, 1, bond);
    bond_graph.add_edge(0, 2, bond);

    Topology topology(bond_graph);
    REQUIRE(topology.num_bonds() == 2);
    REQUIRE(topology.num_angles() >= 1); // Should have at least H-C-H angle
  }
}

TEST_CASE("Bond Management", "[topology][bonds]") {
  // TODO: Catch warnings from spdlog
  Topology topology;

  SECTION("Adding bonds") {
    topology.add_bond(0, 1, 1.5);
    REQUIRE(topology.num_bonds() == 1);
    REQUIRE(topology.has_bond(0, 1));
    REQUIRE(topology.has_bond(1, 0)); // Should work both ways

    // Adding duplicate bond should warn but not add
    topology.add_bond(0, 1, 1.5);
    REQUIRE(topology.num_bonds() == 1);
  }

  SECTION("Removing bonds") {
    topology.add_bond(0, 1, 1.5);
    topology.add_bond(1, 2, 1.5);
    REQUIRE(topology.num_bonds() == 2);

    topology.remove_bond(0, 1);
    REQUIRE(topology.num_bonds() == 1);
    REQUIRE_FALSE(topology.has_bond(0, 1));
    REQUIRE(topology.has_bond(1, 2));

    // Removing non-existent bond should warn but not crash
    topology.remove_bond(0, 1);
    REQUIRE(topology.num_bonds() == 1);
  }

  SECTION("Clearing bonds") {
    topology.add_bond(0, 1, 1.5);
    topology.add_bond(1, 2, 1.5);
    topology.add_bond(2, 3, 1.5);

    topology.clear_bonds();
    REQUIRE(topology.num_bonds() == 0);
    REQUIRE(topology.num_angles() == 0);
    REQUIRE(topology.num_dihedrals() == 0);
  }

  SECTION("Getting bonds") {
    topology.add_bond(0, 1, 1.5);
    topology.add_bond(1, 2, 2.0);

    auto bonds = topology.get_bonds();
    REQUIRE(bonds.size() == 2);

    // Check bond lengths
    bool found_15 = false, found_20 = false;
    for (const auto &bond : bonds) {
      if (Catch::Approx(bond.bond_length) == 1.5)
        found_15 = true;
      if (Catch::Approx(bond.bond_length) == 2.0)
        found_20 = true;
    }
    REQUIRE(found_15);
    REQUIRE(found_20);
  }
}

TEST_CASE("Angle Management", "[topology][angles]") {
  Topology topology;

  SECTION("Adding angles manually") {
    topology.add_angle(0, 1, 2);
    REQUIRE(topology.num_angles() == 1);
    REQUIRE(topology.has_angle(0, 1, 2));

    // Adding duplicate should warn but not add
    topology.add_angle(0, 1, 2);
    REQUIRE(topology.num_angles() == 1);
  }

  SECTION("Removing angles") {
    topology.add_angle(0, 1, 2);
    topology.add_angle(1, 2, 3);
    REQUIRE(topology.num_angles() == 2);

    topology.remove_angle(0, 1, 2);
    REQUIRE(topology.num_angles() == 1);
    REQUIRE_FALSE(topology.has_angle(0, 1, 2));
    REQUIRE(topology.has_angle(1, 2, 3));
  }

  SECTION("Clearing angles") {
    topology.add_angle(0, 1, 2);
    topology.add_angle(1, 2, 3);

    topology.clear_angles();
    REQUIRE(topology.num_angles() == 0);
  }

  SECTION("Getting angles") {
    topology.add_angle(0, 1, 2);
    topology.add_angle(1, 2, 3);

    auto angles = topology.get_angles();
    REQUIRE(angles.size() == 2);

    // Check angle structure
    for (const auto &angle : angles) {
      REQUIRE(angle.atom_indices.size() == 3);
      REQUIRE((angle.center_atom() == 1 || angle.center_atom() == 2));
    }
  }
}

TEST_CASE("Dihedral Management", "[topology][dihedrals]") {
  Topology topology;

  SECTION("Adding dihedrals manually") {
    topology.add_dihedral(0, 1, 2, 3, DihedralType::PROPER);
    REQUIRE(topology.num_dihedrals() == 1);
    REQUIRE(topology.num_proper_dihedrals() == 1);
    REQUIRE(topology.num_improper_dihedrals() == 0);
    REQUIRE(topology.has_dihedral(0, 1, 2, 3));

    topology.add_dihedral(4, 5, 6, 7, DihedralType::IMPROPER);
    REQUIRE(topology.num_dihedrals() == 2);
    REQUIRE(topology.num_proper_dihedrals() == 1);
    REQUIRE(topology.num_improper_dihedrals() == 1);
  }

  SECTION("Removing dihedrals") {
    topology.add_dihedral(0, 1, 2, 3, DihedralType::PROPER);
    topology.add_dihedral(1, 2, 3, 4, DihedralType::PROPER);
    REQUIRE(topology.num_dihedrals() == 2);

    topology.remove_dihedral(0, 1, 2, 3);
    REQUIRE(topology.num_dihedrals() == 1);
    REQUIRE_FALSE(topology.has_dihedral(0, 1, 2, 3));
    REQUIRE(topology.has_dihedral(1, 2, 3, 4));
  }

  SECTION("Clearing dihedrals") {
    topology.add_dihedral(0, 1, 2, 3, DihedralType::PROPER);
    topology.add_dihedral(1, 2, 3, 4, DihedralType::IMPROPER);

    topology.clear_dihedrals();
    REQUIRE(topology.num_dihedrals() == 0);
    REQUIRE(topology.num_proper_dihedrals() == 0);
    REQUIRE(topology.num_improper_dihedrals() == 0);
  }

  SECTION("Getting dihedrals") {
    topology.add_dihedral(0, 1, 2, 3, DihedralType::PROPER);
    topology.add_dihedral(4, 5, 6, 7, DihedralType::IMPROPER);

    auto dihedrals = topology.get_dihedrals();
    REQUIRE(dihedrals.size() == 2);

    bool found_proper = false, found_improper = false;
    for (const auto &dihedral : dihedrals) {
      REQUIRE(dihedral.atom_indices.size() == 4);
      if (dihedral.type == DihedralType::PROPER)
        found_proper = true;
      if (dihedral.type == DihedralType::IMPROPER)
        found_improper = true;
    }
    REQUIRE(found_proper);
    REQUIRE(found_improper);
  }
}

TEST_CASE("Automatic Generation from Bonds", "[topology][generation]") {
  Topology topology;

  SECTION("Generate angles from bonds") {
    // Create a bent molecule: 0-1-2
    topology.add_bond(0, 1, 1.0);
    topology.add_bond(1, 2, 1.0);

    topology.generate_angles_from_bonds();
    REQUIRE(topology.num_angles() == 1);
    REQUIRE(topology.has_angle(0, 1, 2));
  }

  SECTION("Generate proper dihedrals from bonds") {
    // Create a chain: 0-1-2-3
    topology.add_bond(0, 1, 1.0);
    topology.add_bond(1, 2, 1.0);
    topology.add_bond(2, 3, 1.0);

    topology.generate_proper_dihedrals_from_bonds();
    auto dihedrals = topology.get_dihedrals();
    REQUIRE(topology.num_proper_dihedrals() == 1);
    REQUIRE(topology.has_dihedral(0, 1, 2, 3));
  }

  SECTION("Generate improper dihedrals from bonds") {
    // Create a tetrahedral center: 1-0, 2-0, 3-0 (center at 0)
    topology.add_bond(0, 1, 1.0);
    topology.add_bond(0, 2, 1.0);
    topology.add_bond(0, 3, 1.0);

    topology.generate_improper_dihedrals_from_bonds();
    REQUIRE(topology.num_improper_dihedrals() == 1);
  }

  SECTION("Generate all from bonds") {
    auto atoms = create_linear_molecule();
    Topology topology(atoms);

    // Linear H-C-C-H should have:
    // - 3 bonds (H-C, C-C, C-H)
    // - 2 angles (H-C-C, C-C-H)
    // - 1 proper dihedral (H-C-C-H)
    // - 0 improper dihedrals

    REQUIRE(topology.num_bonds() >= 3);
    REQUIRE(topology.num_angles() >= 2);
    REQUIRE(topology.num_proper_dihedrals() >= 1);
  }
}

TEST_CASE("Graph Queries", "[topology][graph]") {
  Topology topology;

  SECTION("Get bonded atoms") {
    topology.add_bond(0, 1, 1.0);
    topology.add_bond(0, 2, 1.0);
    topology.add_bond(0, 3, 1.0);

    auto bonded = topology.get_bonded_atoms(0);
    REQUIRE(bonded.size() == 3);

    std::sort(bonded.begin(), bonded.end());
    REQUIRE(bonded == std::vector<size_t>{1, 2, 3});

    // Non-bonded atom should return empty
    auto empty = topology.get_bonded_atoms(4);
    REQUIRE(empty.empty());
  }

  SECTION("Get atoms at distance") {
    // Create chain: 0-1-2-3-4
    topology.add_bond(0, 1, 1.0);
    topology.add_bond(1, 2, 1.0);
    topology.add_bond(2, 3, 1.0);
    topology.add_bond(3, 4, 1.0);

    // Distance 0
    auto dist0 = topology.get_atoms_at_distance(2, 0);
    REQUIRE(dist0.size() == 1);
    REQUIRE(dist0[0] == 2);

    // Distance 1
    auto dist1 = topology.get_atoms_at_distance(2, 1);
    REQUIRE(dist1.size() == 2);
    std::sort(dist1.begin(), dist1.end());
    REQUIRE(dist1 == std::vector<size_t>{1, 3});

    // Distance 2
    auto dist2 = topology.get_atoms_at_distance(2, 2);
    REQUIRE(dist2.size() == 2);
    std::sort(dist2.begin(), dist2.end());
    REQUIRE(dist2 == std::vector<size_t>{0, 4});

    // Distance 3 (beyond end of chain)
    auto dist3 = topology.get_atoms_at_distance(2, 3);
    REQUIRE(dist3.empty());
  }
}

TEST_CASE("Connected Components and Molecules", "[topology][molecules]") {
  SECTION("Single connected component") {
    auto atoms = create_test_atoms();
    Topology topology(atoms);

    auto molecules = topology.extract_molecules();
    REQUIRE(molecules.size() == 1);
    REQUIRE(molecules[0].size() == 5); // C + 4H
  }

  SECTION("Multiple disconnected components") {
    auto atoms = create_disconnected_atoms();
    Topology topology(atoms);

    auto molecules = topology.extract_molecules();
    REQUIRE(molecules.size() == 2);
    REQUIRE(molecules[0].size() == 2); // H2
    REQUIRE(molecules[1].size() == 2); // H2
  }

  SECTION("Empty topology") {
    std::vector<Atom> empty_atoms;
    Topology topology(empty_atoms);

    auto molecules = topology.extract_molecules();
    REQUIRE(molecules.empty());
  }
}

TEST_CASE("Topology Validation", "[topology][validation]") {
  Topology topology;

  SECTION("Valid topology") {
    topology.add_bond(0, 1, 1.0);
    topology.add_bond(1, 2, 1.0);
    topology.add_angle(0, 1, 2);
    topology.add_dihedral(0, 1, 2, 3, DihedralType::PROPER);

    REQUIRE(topology.validate_topology());
    auto issues = topology.check_issues();
    REQUIRE(issues.empty());
  }

  SECTION("Invalid angles") {
    // Same atom indices in angle
    topology.add_angle(1, 1, 2); // center same as end

    REQUIRE_FALSE(topology.validate_topology());
    auto issues = topology.check_issues();
    REQUIRE_FALSE(issues.empty());
  }

  SECTION("Invalid dihedrals") {
    // Duplicate atom indices in dihedral
    topology.add_dihedral(1, 1, 2, 3, DihedralType::PROPER);

    REQUIRE_FALSE(topology.validate_topology());
    auto issues = topology.check_issues();
    REQUIRE_FALSE(issues.empty());
  }
}

TEST_CASE("Topology String Representation", "[topology][string]") {
  Topology topology;
  topology.add_bond(0, 1, 1.0);
  topology.add_bond(1, 2, 1.0);
  topology.generate_all_from_bonds();

  SECTION("String representation") {
    std::string str = topology.to_string();
    REQUIRE_FALSE(str.empty());
    REQUIRE(str.find("Topology") != std::string::npos);
    REQUIRE(str.find("bonds=") != std::string::npos);
    REQUIRE(str.find("angles=") != std::string::npos);
  }
}

TEST_CASE("Bond Graph Access", "[topology][bondgraph]") {
  auto atoms = create_test_atoms();
  Topology topology(atoms);

  SECTION("Const access to bond graph") {
    const auto &bond_graph = topology.get_bond_graph();
    const auto &adj_list = bond_graph.get_adjacency_list();
    REQUIRE_FALSE(adj_list.empty());
  }

  SECTION("Mutable access to bond graph") {
    auto &bond_graph = topology.get_bond_graph();
    Bond new_bond(2.0);
    bond_graph.add_edge(10, 11, new_bond);

    // Should be able to modify through reference
    const auto &adj_list = bond_graph.get_adjacency_list();
    REQUIRE(adj_list.find(10) != adj_list.end());
  }
}

TEST_CASE("Edge Cases and Error Handling", "[topology][edge_cases]") {
  SECTION("Operations on empty topology") {
    Topology topology;

    REQUIRE(topology.get_bonded_atoms(0).empty());
    REQUIRE(topology.get_atoms_at_distance(0, 1).empty());
    REQUIRE_FALSE(topology.has_bond(0, 1));
    REQUIRE_FALSE(topology.has_angle(0, 1, 2));
    REQUIRE_FALSE(topology.has_dihedral(0, 1, 2, 3));
  }

  SECTION("Large atom indices") {
    Topology topology;
    topology.add_bond(1000, 2000, 1.0);

    REQUIRE(topology.has_bond(1000, 2000));
    REQUIRE(topology.num_bonds() == 1);
  }

  SECTION("Self-loops and invalid operations") {
    Topology topology;

    // Self-bond (though this should be handled by validation)
    topology.add_bond(0, 0, 1.0);
    REQUIRE(topology.num_bonds() == 1);

    // Self-angle (invalid)
    topology.add_angle(0, 0, 1);
    auto issues = topology.check_issues();
    REQUIRE_FALSE(issues.empty());
  }
}

TEST_CASE("Trajectory and Topology Integration", "[trajectory][topology]") {
  SECTION("Load PDB and check topology") {
    Trajectory trajectory;
    std::vector<fs::path> files = {EXAMPLES_DIR / "coord.pdb"};
    trajectory.load_files(files);

    REQUIRE(trajectory.next_frame());

    const Topology &topology = trajectory.get_topology();

    // For a typical water box, we expect many bonds, angles, and molecules
    // The exact numbers depend on the PDB content and bond detection logic
    // Let's assume some reasonable minimums for a non-empty PDB
    REQUIRE(topology.num_bonds() > 0);
    REQUIRE(topology.num_angles() > 0);
    // Dihedrals might be 0 for simple water, but good to check if any are found
    REQUIRE(topology.num_dihedrals() >= 0);

    auto molecules = topology.extract_molecules();
    REQUIRE(molecules.size() > 0); // Should find water molecules

    // Check if the number of atoms in the topology matches the frame
    REQUIRE(topology.num_atoms() == trajectory.num_atoms());

    auto molecules2 = trajectory.extract_molecules();
    REQUIRE(molecules2.size() > 0); // Should find water molecules
    REQUIRE(molecules2.size() ==
            molecules.size()); // Should find water molecules
  }
}

TEST_CASE("Performance and Memory", "[topology][performance]") {
  SECTION("Large number of bonds") {
    Topology topology;

    // Add many bonds
    const size_t num_bonds = 1000;
    for (size_t i = 0; i < num_bonds; ++i) {
      topology.add_bond(i, i + num_bonds, 1.0);
    }

    REQUIRE(topology.num_bonds() == num_bonds);

    // Clear should work efficiently
    topology.clear_bonds();
    REQUIRE(topology.num_bonds() == 0);
  }

  SECTION("Complex molecule generation") {
    // Create a more complex branched molecule
    std::vector<Atom> atoms;
    for (int i = 0; i < 20; ++i) {
      atoms.emplace_back(Vec3{static_cast<double>(i), 0.0, 0.0}, "C", i);
    }

    Topology topology(atoms);
    topology.generate_all_from_bonds();

    // Should handle generation without issues
    REQUIRE(topology.validate_topology());
  }
}
