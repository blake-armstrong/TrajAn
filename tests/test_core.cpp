#include <catch2/benchmark/catch_benchmark.hpp>
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
#include <trajan/core/trajectory.h>
#include <trajan/core/unit_cell.h>
#include <trajan/core/units.h>
#include <trajan/io/selection.h>
#include <vector>

using trajan::Mat3N;
using trajan::Vec3;
using trajan::core::Atom;
using trajan::core::Cell;
using trajan::core::CellList;
using trajan::core::Entity;
using trajan::core::Frame;
using trajan::core::Molecule;
using trajan::core::Trajectory;
using trajan::core::UnitCell;
using trajan::core::VariantEqual;
using trajan::core::VariantHash;
using trajan::core::element::Element;
using trajan::io::SelectionCriteria;
using trajan::io::SelectionParser;
using Entities = std::vector<trajan::core::EntityType>;

namespace units = trajan::units;

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
  auto result = trajan::util::combine_deduplicate_map(
      all_entities, VariantHash(), VariantEqual());
  Entities deduplicated_entities = result.first;
  std::vector<std::bitset<8>> presence_tracker = result.second;
  size_t entities_size = deduplicated_entities.size();

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
  SECTION("Comparing pair counting between methods with presence tracker") {
    std::atomic<size_t> cell_list_pairs{0};
    std::atomic<size_t> cell_list_pairs_presence{0};

    trajan::core::NeighbourList neighbour_list(unit_cell, cutoff, num_threads);
    neighbour_list.update(deduplicated_entities);
    neighbour_list.iterate_neighbours([&](const trajan::core::Entity &e1,
                                          const trajan::core::Entity &e2,
                                          double rsq) { cell_list_pairs++; });
    neighbour_list.iterate_neighbours([&](const trajan::core::Entity &e1,
                                          const trajan::core::Entity &e2,
                                          double rsq) {
      std::bitset<8> source1 = presence_tracker[e1.idx];
      std::bitset<8> source2 = presence_tracker[e2.idx];
      bool is_cross_section = (source1.test(0) && source2.test(1)) ||
                              (source1.test(1) && source2.test(0));
      if (!is_cross_section) {
        return;
      }
      cell_list_pairs_presence++;
    });
    REQUIRE(cell_list_pairs == cell_list_pairs_presence);
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

// TEST_CASE("Molecule identification from atoms", "[molecule]") {
//   SECTION("Empty atom list") {
//     std::vector<Atom> atoms;
//     auto molecules = trajan::core::identify_molecules(atoms);
//
//     REQUIRE(molecules.empty());
//   }
//
//   SECTION("Single water molecule") {
//     std::vector<Atom> atoms = {
//         Atom(Vec3{0.0, 0.0, 0.0}, Element("O", true), 0),
//         Atom(Vec3{0.8, 0.0, 0.0}, Element("H", true), 1),
//         Atom(Vec3{-0.2, 0.8, 0.0}, Element("H", true), 2)};
//
//     auto molecules = trajan::core::identify_molecules(atoms);
//
//     REQUIRE(molecules.size() == 1);
//     REQUIRE(molecules[0].size() == 3);
//   }
//
//   SECTION("Two separate water molecules") {
//     std::vector<Atom> atoms = {
//         // First water molecule
//         Atom(Vec3{0.0, 0.0, 0.0}, Element("O", true), 0),
//         Atom(Vec3{0.8, 0.0, 0.0}, Element("H", true), 1),
//         Atom(Vec3{-0.2, 0.8, 0.0}, Element("H", true), 2),
//         // Second water molecule (far away)
//         Atom(Vec3{10.0, 10.0, 10.0}, Element("O", true), 4),
//         Atom(Vec3{10.8, 10.0, 10.0}, Element("H", true), 5),
//         Atom(Vec3{9.8, 10.8, 10.0}, Element("H", true), 6)};
//
//     auto molecules = trajan::core::identify_molecules(atoms);
//
//     REQUIRE(molecules.size() == 2);
//     REQUIRE(molecules[0].size() == 3);
//     REQUIRE(molecules[1].size() == 3);
//   }
//
//   SECTION("Three separate atoms that are not bonded") {
//     std::vector<Atom> atoms = {
//         Atom(Vec3{0.0, 0.0, 0.0}, Element("Na", true), 0),
//         Atom(Vec3{10.0, 0.0, 0.0}, Element("Cl", true), 1),
//         Atom(Vec3{0.0, 10.0, 0.0}, Element("Ar", true), 2)};
//
//     auto molecules = trajan::core::identify_molecules(atoms);
//
//     REQUIRE(molecules.size() == 3);
//     REQUIRE(molecules[0].size() == 1);
//     REQUIRE(molecules[1].size() == 1);
//     REQUIRE(molecules[2].size() == 1);
//   }
//
//   SECTION("Custom bond tolerance") {
//     std::vector<Atom> atoms = {
//         Atom(Vec3{0.0, 0.0, 0.0}, Element("C", true), 0),
//         Atom(Vec3{2.0, 0.0, 0.0}, Element("C", true), 1)
//         // Distance is larger than default tolerance
//     };
//
//     // With default tolerance, should be two separate molecules
//     auto molecules_default = trajan::core::identify_molecules(atoms);
//     REQUIRE(molecules_default.size() == 2);
//
//     // With increased tolerance, should be one molecule
//     auto molecules_increased = trajan::core::identify_molecules(atoms, 1.0);
//     REQUIRE(molecules_increased.size() == 1);
//     REQUIRE(molecules_increased[0].size() == 2);
//   }
// }

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

// TEST_CASE("Complex molecule identification", "[molecule][integration]") {
//   SECTION("Methane molecule") {
//     std::vector<Atom> atoms = {
//         Atom(Vec3{0.0, 0.0, 0.0}, Element("C", true), 0),
//         Atom(Vec3{0.8, 0.8, 0.8}, Element("H", true), 1),
//         Atom(Vec3{-0.8, -0.8, 0.8}, Element("H", true), 2),
//         Atom(Vec3{0.8, -0.8, -0.8}, Element("H", true), 3),
//         Atom(Vec3{-0.8, 0.8, -0.8}, Element("H", true), 4)};
//
//     auto molecules = trajan::core::identify_molecules(atoms, 0.8);
//
//     REQUIRE(molecules.size() == 1);
//     REQUIRE(molecules[0].size() == 5);
//   }
//
//   SECTION("Ethanol molecule") {
//     std::vector<Atom> atoms = {
//         // Carbon chain
//         Atom(Vec3{0.0, 0.0, 0.0}, Element("C", true), 0),
//         Atom(Vec3{1.5, 0.0, 0.0}, Element("C", true), 1),
//         // Hydrogens on first carbon
//         Atom(Vec3{-0.5, 0.9, 0.0}, Element("H", true), 2),
//         Atom(Vec3{-0.5, -0.5, 0.9}, Element("H", true), 3),
//         Atom(Vec3{-0.5, -0.5, -0.9}, Element("H", true), 4),
//         // Hydrogens on second carbon
//         Atom(Vec3{2.0, 0.9, 0.0}, Element("H", true), 5),
//         Atom(Vec3{2.0, -0.5, -0.9}, Element("H", true), 6),
//         // Oxygen and its hydrogen
//         Atom(Vec3{2.0, -0.5, 0.9}, Element("O", true), 7),
//         Atom(Vec3{3.0, -0.5, 0.9}, Element("H", true), 8)};
//
//     auto molecules = trajan::core::identify_molecules(atoms);
//
//     REQUIRE(molecules.size() == 1);
//     REQUIRE(molecules[0].size() == 9);
//   }
//
//   SECTION("Mixture of molecules") {
//     std::vector<Atom> atoms = {
//         // Water molecule
//         Atom(Vec3{0.0, 0.0, 0.0}, Element("O", true), 0),
//         Atom(Vec3{0.8, 0.0, 0.0}, Element("H", true), 1),
//         Atom(Vec3{-0.2, 0.8, 0.0}, Element("H", true), 2),
//         // CO2 molecule
//         Atom(Vec3{10.0, 10.0, 10.0}, Element("C", true), 3),
//         Atom(Vec3{11.2, 10.0, 10.0}, Element("O", true), 4),
//         Atom(Vec3{8.8, 10.0, 10.0}, Element("O", true), 5),
//         // Single atom
//         Atom(Vec3{20.0, 20.0, 20.0}, Element("Na", true), 6)};
//
//     auto molecules = trajan::core::identify_molecules(atoms);
//
//     REQUIRE(molecules.size() == 3);
//     // Check sizes of each molecule
//     bool found_water = false;
//     bool found_co2 = false;
//     bool found_sodium = false;
//
//     for (const auto &mol : molecules) {
//       if (mol.size() == 3) {
//         // Could be water or CO2
//         // In a real test, you'd need more checks to distinguish
//         if (!found_water) {
//           found_water = true;
//         } else {
//           found_co2 = true;
//         }
//       } else if (mol.size() == 1) {
//         found_sodium = true;
//       }
//     }
//
//     REQUIRE(found_water);
//     REQUIRE(found_co2);
//     REQUIRE(found_sodium);
//   }
// }
