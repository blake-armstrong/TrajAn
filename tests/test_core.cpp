// tests
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <chrono>
#include <iostream>
#include <omp.h>
#include <random>
#include <trajan/core/element.h>
#include <trajan/core/linear_algebra.h>
#include <trajan/core/neigh.h>
#include <trajan/core/unit_cell.h>
#include <trajan/core/units.h>
#include <vector>

using trajan::Mat3N;
using trajan::Vec3;
using trajan::core::Atom;
using trajan::core::Cell;
using trajan::core::CellList;
using trajan::core::Element;
using trajan::core::UnitCell;

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

// Helper function to generate random atoms within a unit cell
// std::vector<Atom> generate_random_atoms(const UnitCell &cell, int num_atoms)
// {
//   std::vector<Atom> atoms;
//   std::random_device rd;
//   std::mt19937 gen(rd());
//   std::uniform_real_distribution<> dis(0.0, 1.0);
//
//   for (int i = 0; i < num_atoms; ++i) {
//     Vec3 position(dis(gen), dis(gen), dis(gen));
//     position = cell.to_cartesian(position);
//     atoms.emplace_back(position, 0, i);
//   }
//   return atoms;
// }

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
    cell.add_atom(atom);
    REQUIRE(cell.get_atoms().size() == 1);
    REQUIRE(cell.get_atoms()[0].position() == pos);
  }

  SECTION("Adding multiple atoms") {
    const int num_atoms = 10;
    std::vector<Atom> atoms;
    for (int i = 0; i < num_atoms; ++i) {
      atoms.emplace_back(Vec3(i, i, i), 0, i);
      cell.add_atom(atoms.back());
    }
    REQUIRE(cell.get_atoms().size() == num_atoms);
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
  const double box_size = 25.0;
  UnitCell unit_cell = trajan::core::cubic_cell(box_size);
  const int num_atoms = 10000;
  const double cutoff = 9.0;
  const int num_threads = 1;
  std::vector<Atom> atoms;
  atoms.reserve(num_atoms);
  std::random_device rd;
  std::mt19937 gen(rd());
  // std::mt19937 gen(42);
  std::uniform_real_distribution<> dis(-5, box_size + 5);

  // create atoms with random positions
  Mat3N atoms_pos(3, num_atoms);
  for (int i = 0; i < num_atoms; ++i) {
    Vec3 position(dis(gen), dis(gen), dis(gen));
    atoms.emplace_back(position, 0, i);
    atoms_pos(0, i) = position.x();
    atoms_pos(1, i) = position.y();
    atoms_pos(2, i) = position.z();
  }

  SECTION("Comparing pair counting between methods") {
    std::atomic<size_t> cell_list_pairs{0};
    std::atomic<size_t> double_loop_pairs{0};

    // count pairs using cell list
    CellList cell_list(unit_cell, cutoff, num_threads);
    double cell_list_time = measure_execution_time([&]() {
      cell_list.update(atoms, atoms_pos);
      cell_list.for_each_pair([&](const Atom &a1, const Atom &a2, double rsq) {
        cell_list_pairs++;
      });
    });

    // count pairs using double loop
    double double_loop_time = measure_execution_time([&]() {
      for (size_t i = 0; i < atoms.size(); ++i) {
        for (size_t j = i + 1; j < atoms.size(); ++j) {
          Vec3 dr = atoms[j].position() - atoms[i].position();
          Vec3 frac_dr = unit_cell.to_fractional(dr);
          frac_dr = frac_dr.array() - frac_dr.array().round();
          dr = unit_cell.to_cartesian(frac_dr);

          if (dr.norm() <= cutoff) {
            double_loop_pairs++;
          }
        }
      }
    });

    INFO("System size: " << num_atoms << " atoms");
    INFO("Box size: " << box_size << "x" << box_size << "x" << box_size);
    INFO("Cutoff distance: " << cutoff);
    INFO("Number of threads: " << num_threads);
    // estimate the expected number of pairs
    double volume = box_size * box_size * box_size;
    double expected_pairs_per_atom =
        (4.0 / 3.0) * units::PI * std::pow(cutoff, 3) / volume * num_atoms;
    double total_expected_pairs = expected_pairs_per_atom * num_atoms / 2;
    INFO("Expected approximate number of pairs: ~" << total_expected_pairs);
    INFO("Pairs found (cell list): " << cell_list_pairs);
    INFO("Pairs found (double loop): " << double_loop_pairs);
    // Verify that both methods find the same number of pairs
    REQUIRE(cell_list_pairs == double_loop_pairs);
    REQUIRE(static_cast<double>(cell_list_pairs) > total_expected_pairs * 0.8);
    REQUIRE(static_cast<double>(cell_list_pairs) < total_expected_pairs * 1.2);
    INFO("Cell list execution time: " << cell_list_time << " seconds");
    INFO("Double loop execution time: " << double_loop_time << " seconds");
    INFO("Speed-up factor: " << double_loop_time / cell_list_time);
    // cell list should be faster
    REQUIRE(double_loop_time > cell_list_time);
  }
}

// TEST_CASE("CellList Performance", "[cell_list_performance]") {
//   const double box_size = 25.0;
//   const double cutoff = 9.0;
//   const int num_atoms = 5000;
//
//   // Create unit cell and atoms
//   UnitCell unit_cell = trajan::core::cubic_cell(box_size);
//   std::vector<Atom> atoms;
//   atoms.reserve(num_atoms);
//
//   // Create deterministic random number generator for reproducible tests
//   std::mt19937 gen(42);
//   std::uniform_real_distribution<> dis(0.0, box_size);
//
//   for (int i = 0; i < num_atoms; ++i) {
//     Vec3 position(dis(gen), dis(gen), dis(gen));
//     atoms.emplace_back(position, 0, i);
//   }
//
//   SECTION("Correctness Test") {
//     // Create two cell lists with different thread counts
//     CellList cell_list_single(unit_cell, cutoff, 1);
//     CellList cell_list_multi(unit_cell, cutoff, 4);
//
//     // Update both cell lists
//     cell_list_single.update(atoms);
//     cell_list_multi.update(atoms);
//
//     // Count pairs in both cell lists
//     std::atomic<size_t> pairs_single{0};
//     std::atomic<size_t> pairs_multi{0};
//     std::vector<std::pair<int, int>> pair_indices_single;
//     std::vector<std::pair<int, int>> pair_indices_multi;
//
//     cell_list_single.for_each_pair([&](const Atom &a1, const Atom &a2) {
//       pairs_single++;
//       pair_indices_single.emplace_back(a1.id(), a2.id());
//     });
//
//     cell_list_multi.for_each_pair([&](const Atom &a1, const Atom &a2) {
//       pairs_multi++;
//       pair_indices_multi.emplace_back(a1.id(), a2.id());
//     });
//
//     // Sort pair indices for comparison
//     std::sort(pair_indices_single.begin(), pair_indices_single.end());
//     std::sort(pair_indices_multi.begin(), pair_indices_multi.end());
//
//     // Verify results are identical
//     REQUIRE(pairs_single == pairs_multi);
//     REQUIRE(pair_indices_single == pair_indices_multi);
//   }
//
//   // SECTION("Performance Test") {
//   //   const int num_runs = 5;
//   //   std::vector<double> single_thread_times;
//   //   std::vector<double> multi_thread_times;
//   //
//   //   CellList cell_list_single(unit_cell, cutoff, 1);
//   //   CellList cell_list_multi(unit_cell, cutoff, 4);
//   //
//   //   for (int i = 0; i < num_runs; ++i) {
//   //     // Measure single-threaded performance
//   //     auto start_single = std::chrono::high_resolution_clock::now();
//   //     cell_list_single.update(atoms);
//   //     auto end_single = std::chrono::high_resolution_clock::now();
//   //     single_thread_times.push_back(
//   //         std::chrono::duration<double>(end_single -
//   start_single).count());
//   //
//   //     // Measure multi-threaded performance
//   //     auto start_multi = std::chrono::high_resolution_clock::now();
//   //     cell_list_multi.update(atoms);
//   //     auto end_multi = std::chrono::high_resolution_clock::now();
//   //     multi_thread_times.push_back(
//   //         std::chrono::duration<double>(end_multi - start_multi).count());
//   //   }
//   //
//   //   // Calculate average times
//   //   double avg_single = std::accumulate(single_thread_times.begin(),
//   //                                       single_thread_times.end(), 0.0) /
//   //                       num_runs;
//   //   double avg_multi = std::accumulate(multi_thread_times.begin(),
//   //                                      multi_thread_times.end(), 0.0) /
//   //                      num_runs;
//   //
//   //   // Verify multi-threaded version is faster
//   //   INFO("Average single-threaded time: " << avg_single << " seconds");
//   //   INFO("Average multi-threaded time: " << avg_multi << " seconds");
//   //   INFO("Speedup factor: " << avg_single / avg_multi);
//   //
//   //   REQUIRE(avg_multi < avg_single);
//   // }
// }
