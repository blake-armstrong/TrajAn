#pragma once
#include <CLI/CLI.hpp>
#include <Eigen/src/Geometry/AngleAxis.h>
#include <Eigen/src/Geometry/Quaternion.h>
#include <nlohmann/json.hpp>
#include <occ/core/linear_algebra.h>
#include <trajan/core/topology.h>
#include <trajan/core/trajectory.h>

namespace trajan::main {

using Eigen::AngleAxisd;
using Eigen::Quaterniond;
constexpr double pi = occ::constants::pi<double>;
using nlohmann::json;
using trajan::core::Atom;
using trajan::core::ATOM_RESTRICTIONS;
using trajan::core::Trajectory;

struct QCOpts {
  std::string carbonate_selection;
  bool use_c3_symmetry{true};
  size_t n_clusters{3};
  std::string output_file{"qc.out"};
  bool train_mode{false};
  std::string reference_config_file{"surface_configs.json"};
  double distance_threshold{0.3};
  int max_kmeans_iter{100};
};

struct SurfaceState {
  std::vector<Quaterniond> carbonate_orientations;
  size_t frame_number;

  inline occ::DVec to_feature_vector() const {
    int n = carbonate_orientations.size();
    occ::DVec features(4 * n);

    for (int i = 0; i < n; ++i) {
      const auto &q = carbonate_orientations[i];
      features(4 * i + 0) = q.w();
      features(4 * i + 1) = q.x();
      features(4 * i + 2) = q.y();
      features(4 * i + 3) = q.z();
    }

    return features;
  }
};

json serialize_surface_state(const SurfaceState &state);
SurfaceState deserialize_surface_state(const json &j);
std::vector<SurfaceState> load_reference_configs(const std::string &filename);
void save_reference_configs(const std::vector<SurfaceState> &configs,
                            const std::string &filename);

Quaterniond
compute_carbonate_orientation(const core::Molecule &carb,
                              const std::optional<occ::crystal::UnitCell> &uc);
Quaterniond apply_c3_symmetry(const Quaterniond &q);
double surface_state_distance(const SurfaceState &s1, const SurfaceState &s2);
double surface_state_distance_with_permutation(const SurfaceState &s1,
                                               const SurfaceState &s2);
std::vector<SurfaceState>
kmeans_surface_states(const std::vector<SurfaceState> &states, int k,
                      int max_iter = 100);

void run_qc_train(QCOpts const &opts, Trajectory &traj);
void run_qc_analyse(QCOpts const &opts, Trajectory &traj);
void run_qc_subcommand(QCOpts const &opts, Trajectory &traj);
CLI::App *add_qc_subcommand(CLI::App &app, Trajectory &traj);

} // namespace trajan::main
