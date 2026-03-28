#include <stdexcept>
#include <trajan/core/neigh.h>
#include <trajan/core/topology.h>
#include <trajan/core/util.h>
#include <trajan/io/selection.h>
#include <trajan/main/trajan_qc.h>

namespace trajan::main {

json serialize_surface_state(const SurfaceState &state) {
  json j;
  j["frame_number"] = state.frame_number;
  j["n_carbonates"] = state.carbonate_orientations.size();

  json quaternions = json::array();
  for (const auto &q : state.carbonate_orientations) {
    quaternions.push_back(
        {{"w", q.w()}, {"x", q.x()}, {"y", q.y()}, {"z", q.z()}});
  }
  j["quaternions"] = quaternions;

  return j;
}

// Deserialize from JSON
SurfaceState deserialize_surface_state(const json &j) {
  SurfaceState state;
  state.frame_number = j["frame_number"];

  for (const auto &q_json : j["quaternions"]) {
    Quaterniond q(q_json["w"], q_json["x"], q_json["y"], q_json["z"]);
    state.carbonate_orientations.push_back(q);
  }

  return state;
}

// Save reference configurations to file
void save_reference_configs(const std::vector<SurfaceState> &configs,
                            const std::string &filename) {
  nlohmann::json j;
  j["n_configs"] = configs.size();
  j["n_carbonates"] = configs[0].carbonate_orientations.size();

  nlohmann::json configs_json = nlohmann::json::array();
  for (size_t i = 0; i < configs.size(); ++i) {
    configs_json.push_back(serialize_surface_state(configs[i]));
  }
  j["configurations"] = configs_json;

  std::ofstream file(filename);
  file << j.dump(2); // Pretty print with indent=2
  file.close();

  trajan::log::info("Saved {} reference configurations to {}", configs.size(),
                    filename);
}

// Load reference configurations from file
std::vector<SurfaceState> load_reference_configs(const std::string &filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Cannot open reference config file: " + filename);
  }

  nlohmann::json j;
  file >> j;

  std::vector<SurfaceState> configs;
  for (const auto &config_json : j["configurations"]) {
    configs.push_back(deserialize_surface_state(config_json));
  }

  trajan::log::info("Loaded {} reference configurations from {}",
                    configs.size(), filename);

  return configs;
}

Quaterniond
compute_carbonate_orientation(const core::Molecule &carb,
                              const std::optional<occ::crystal::UnitCell> &uc) {
  occ::Vec3 C_pos;
  std::vector<occ::Vec3> O_positions;

  for (const Atom &atom : carb.atoms()) {
    auto pos = atom.position();

    if (atom.element.symbol() == "C") {
      C_pos = pos;
    } else if (atom.element.symbol() == "O") {
      O_positions.push_back(pos);
    } else {
      throw std::runtime_error("Non O/C element in molecule.");
    }
  }

  if (O_positions.size() != 3) {
    throw std::runtime_error(
        fmt::format("Carbonate molecule doesn't have 3 oxygens, has {}",
                    O_positions.size()));
  }

  // CRITICAL: Sort oxygens by angle around C to ensure consistent v1, v2
  // ordering This prevents cross product from flipping when PBC wrapping
  // reorders atoms
  std::sort(O_positions.begin(), O_positions.end(),
            [&C_pos, &uc](const occ::Vec3 &a, const occ::Vec3 &b) {
              occ::Vec3 ca = a - C_pos;
              occ::Vec3 cb = b - C_pos;
              if (uc) {
                ca = util::wrap_distance(ca, uc.value());
                cb = util::wrap_distance(cb, uc.value());
              }
              // Sort by angle in xy-plane for consistent ordering
              double angle_a = std::atan2(ca(1), ca(0));
              double angle_b = std::atan2(cb(1), cb(0));
              return angle_a < angle_b;
            });

  // z-axis: normal to carbonate plane
  occ::Vec3 v1 = O_positions[1] - O_positions[0];
  occ::Vec3 v2 = O_positions[2] - O_positions[0];
  if (uc) {
    v1 = util::wrap_distance(v1, uc.value());
    v2 = util::wrap_distance(v2, uc.value());
  }
  occ::Vec3 z_current = v1.cross(v2);

  if (z_current.norm() < 1e-6) {
    trajan::log::warn("Degenerate carbonate geometry");
    return Quaterniond::Identity();
  }

  z_current.normalize();

  // Ensure normal points upward
  if (z_current(2) < 0) {
    z_current = -z_current;
  }

  // x-axis: from C to first oxygen, perpendicular to z
  occ::Vec3 CO = O_positions[0] - C_pos;
  if (uc) {
    CO = util::wrap_distance(CO, uc.value());
  }
  occ::Vec3 x_current = CO - CO.dot(z_current) * z_current;

  if (x_current.norm() < 1e-6) {
    trajan::log::warn("Carbonate has degenerate x-axis");
    return Quaterniond::Identity();
  }

  x_current.normalize();
  occ::Vec3 y_current = z_current.cross(x_current);

  occ::Mat3 R_current;
  R_current.col(0) = x_current;
  R_current.col(1) = y_current;
  R_current.col(2) = z_current;

  return Quaterniond(R_current);
}

// Quaterniond
// compute_carbonate_orientation(const core::Molecule &carb,
//                               const std::optional<occ::crystal::UnitCell>
//                               &uc) {
//   occ::Vec3 C_pos;
//   std::vector<occ::Vec3> O_positions;
//
//   for (const Atom &atom : carb.atoms()) {
//     auto pos = atom.position();
//     trajan::log::info(fmt::format("{} {} {}", pos[0], pos[1], pos[2]));
//     if (atom.element.symbol() == "C") {
//       C_pos = pos;
//     } else if (atom.element.symbol() == "O") {
//       O_positions.push_back(pos);
//     } else {
//       throw std::runtime_error("Non O/C element in molecule.");
//     }
//   }
//
//   if (O_positions.size() != 3) {
//     throw std::runtime_error(
//         fmt::format("Carbonate molecule doesn't have 3 oxygens, has {}",
//                     O_positions.size()));
//   }
//
//   // z-axis: normal to carbonate plane
//   occ::Vec3 v1 = O_positions[1] - O_positions[0];
//   occ::Vec3 v2 = O_positions[2] - O_positions[0];
//   if (uc) {
//     v1 = util::wrap_distance(v1, uc.value());
//     v2 = util::wrap_distance(v2, uc.value());
//   }
//   trajan::log::info("v1: {} {} {}", v1[0], v1[1], v1[2]);
//   trajan::log::info("v2: {} {} {}", v2[0], v2[1], v2[2]);
//   occ::Vec3 z_current = v1.cross(v2);
//
//   if (z_current.norm() < 1e-6) {
//     trajan::log::warn("Degenerate carbonate geometry");
//     return Quaterniond::Identity();
//   }
//
//   z_current.normalize();
//
//   // x-axis: from C to first oxygen, perpendicular to z
//   occ::Vec3 CO = O_positions[0] - C_pos;
//   if (uc) {
//     CO = util::wrap_distance(CO, uc.value());
//   }
//
//   trajan::log::info("CO: {} {} {}", CO[0], CO[1], CO[2]);
//   if (z_current.dot(CO) < 0.0) {
//     z_current = -z_current;
//   }
//   occ::Vec3 x_current = CO - CO.dot(z_current) * z_current;
//
//   if (x_current.norm() < 1e-6) {
//     trajan::log::warn("Carbonate has degenerate x-axis");
//     return Quaterniond::Identity();
//   }
//
//   x_current.normalize();
//   occ::Vec3 y_current = z_current.cross(x_current);
//
//   occ::Mat3 R_current;
//   R_current.col(0) = x_current;
//   R_current.col(1) = y_current;
//   R_current.col(2) = z_current;
//
//   return Quaterniond(R_current);
// }

// Quaterniond apply_c3_symmetry(const Quaterniond &q) {
//   Quaterniond identity(1, 0, 0, 0);
//
//   AngleAxisd rot120(2.0 * pi / 3.0, occ::Vec3::UnitZ());
//   AngleAxisd rot240(4.0 * pi / 3.0, occ::Vec3::UnitZ());
//
//   Quaterniond q_120 = q * Quaterniond(rot120);
//   Quaterniond q_240 = q * Quaterniond(rot240);
//
//   auto distance = [&identity](const Quaterniond &q1) {
//     double dot = std::abs(q1.dot(identity));
//     return 2.0 * std::acos(std::min(1.0, dot));
//   };
//
//   double d0 = distance(q);
//   double d120 = distance(q_120);
//   double d240 = distance(q_240);
//
//   if (d0 <= d120 && d0 <= d240)
//     return q;
//   if (d120 <= d240)
//     return q_120;
//   return q_240;
// }

Quaterniond apply_c3_symmetry(const Quaterniond &q_in) {
  Quaterniond q = q_in.normalized();

  AngleAxisd rot120(2.0 * pi / 3.0, occ::Vec3::UnitZ());
  AngleAxisd rot240(4.0 * pi / 3.0, occ::Vec3::UnitZ());

  Quaterniond q0 = q;
  Quaterniond q120 = q * Quaterniond(rot120);
  Quaterniond q240 = q * Quaterniond(rot240);

  auto score = [](const Quaterniond &qq) { return std::abs(qq.w()); };

  double s0 = score(q0);
  double s120 = score(q120);
  double s240 = score(q240);

  if (s0 >= s120 && s0 >= s240)
    return q0;
  if (s120 >= s240)
    return q120;
  return q240;
}

double surface_state_distance(const SurfaceState &s1, const SurfaceState &s2) {
  if (s1.carbonate_orientations.size() != s2.carbonate_orientations.size()) {
    throw std::runtime_error(
        "Surface states have different number of carbonates");
  }

  double total_dist_sq = 0.0;
  int n = s1.carbonate_orientations.size();

  for (int i = 0; i < n; ++i) {
    double dot = std::abs(
        s1.carbonate_orientations[i].dot(s2.carbonate_orientations[i]));
    dot = std::min(1.0, std::max(-1.0, dot));
    double angle_dist = 2.0 * std::acos(dot);
    total_dist_sq += angle_dist * angle_dist;
  }

  return std::sqrt(total_dist_sq / n); // RMS distance
}

// double surface_state_distance_with_permutation(const SurfaceState &s1,
//                                                const SurfaceState &s2) {
//   int n = s1.carbonate_orientations.size();
//
//   if (n != s2.carbonate_orientations.size()) {
//     throw std::runtime_error(
//         "Surface states have different number of carbonates");
//   }
//
//   // For 4 carbonates, there are 4! = 24 permutations
//   // We'll try all of them and find the minimum distance
//
//   std::vector<int> indices(n);
//   for (int i = 0; i < n; ++i) {
//     indices[i] = i;
//   }
//
//   double min_distance = 1e9;
//
//   // Try all permutations
//   do {
//     double total_dist_sq = 0.0;
//
//     for (int i = 0; i < n; ++i) {
//       int j = indices[i]; // s1[i] matches with s2[j]
//
//       double dot = std::abs(
//           s1.carbonate_orientations[i].dot(s2.carbonate_orientations[j]));
//       dot = std::min(1.0, std::max(-1.0, dot));
//       double angle_dist = 2.0 * std::acos(dot);
//       total_dist_sq += angle_dist * angle_dist;
//     }
//
//     double distance = std::sqrt(total_dist_sq / n);
//     min_distance = std::min(min_distance, distance);
//
//   } while (std::next_permutation(indices.begin(), indices.end()));
//
//   return min_distance;
// }

double surface_state_distance_with_permutation(const SurfaceState &s1,
                                               const SurfaceState &s2) {
  int n = s1.carbonate_orientations.size();

  if (n != s2.carbonate_orientations.size()) {
    throw std::runtime_error(
        "Surface states have different number of carbonates");
  }

  std::vector<int> indices(n);
  for (int i = 0; i < n; ++i) {
    indices[i] = i;
  }

  // Pre-compute C3 rotations for comparison
  AngleAxisd rot120(2.0 * pi / 3.0, occ::Vec3::UnitZ());
  AngleAxisd rot240(4.0 * pi / 3.0, occ::Vec3::UnitZ());
  Quaterniond q_rot120(rot120);
  Quaterniond q_rot240(rot240);

  double min_distance = 1e9;

  // Try all permutations
  do {
    double total_dist_sq = 0.0;

    for (int i = 0; i < n; ++i) {
      int j = indices[i]; // s1[i] matches with s2[j]

      const Quaterniond &q1 = s1.carbonate_orientations[i];
      const Quaterniond &q2 = s2.carbonate_orientations[j];

      // Try all C3-equivalent versions of q2
      Quaterniond q2_0 = q2;
      Quaterniond q2_120 = q2 * q_rot120;
      Quaterniond q2_240 = q2 * q_rot240;

      // Find minimum distance among C3-equivalent quaternions
      double min_angle = 1e9;
      for (const auto &q2_variant : {q2_0, q2_120, q2_240}) {
        double dot = std::abs(q1.dot(q2_variant));
        dot = std::min(1.0, std::max(-1.0, dot));
        double angle = 2.0 * std::acos(dot);
        min_angle = std::min(min_angle, angle);
      }

      total_dist_sq += min_angle * min_angle;
    }

    double distance = std::sqrt(total_dist_sq / n);
    min_distance = std::min(min_distance, distance);

  } while (std::next_permutation(indices.begin(), indices.end()));

  return min_distance;
}

// K-means clustering for surface configurations
std::vector<SurfaceState>
kmeans_surface_states(const std::vector<SurfaceState> &states, int k,
                      int max_iter) {
  int n_states = states.size();
  std::vector<SurfaceState> centroids;

  // Initialize with random states from trajectory
  std::srand(42);
  for (int i = 0; i < k; ++i) {
    int idx = (i * n_states) / k; // Spread initial centroids
    centroids.push_back(states[idx]);
  }

  std::vector<int> assignments(n_states, -1);

  for (int iter = 0; iter < max_iter; ++iter) {
    bool changed = false;

    // Assignment step
    for (int i = 0; i < n_states; ++i) {
      double min_dist = 1e9;
      int best_cluster = 0;

      for (int c = 0; c < k; ++c) {
        double dist = surface_state_distance(states[i], centroids[c]);
        if (dist < min_dist) {
          min_dist = dist;
          best_cluster = c;
        }
      }

      if (assignments[i] != best_cluster) {
        changed = true;
        assignments[i] = best_cluster;
      }
    }

    if (!changed) {
      trajan::log::info("K-means converged after {} iterations", iter);
      break;
    }

    // Update centroids
    for (int c = 0; c < k; ++c) {
      std::vector<SurfaceState> cluster_members;
      for (int i = 0; i < n_states; ++i) {
        if (assignments[i] == c) {
          cluster_members.push_back(states[i]);
        }
      }

      if (cluster_members.empty())
        continue;

      int n_carbs = centroids[c].carbonate_orientations.size();
      SurfaceState new_centroid;
      new_centroid.frame_number = centroids[c].frame_number;

      // Average each carbonate's quaternion
      for (int carb_id = 0; carb_id < n_carbs; ++carb_id) {
        Eigen::Vector4d sum(0, 0, 0, 0);

        for (const auto &state : cluster_members) {
          const auto &q = state.carbonate_orientations[carb_id];
          Eigen::Vector4d q_vec(q.w(), q.x(), q.y(), q.z());

          // Handle quaternion double cover
          if (sum.dot(q_vec) < 0) {
            q_vec = -q_vec;
          }

          sum += q_vec;
        }

        sum.normalize();
        new_centroid.carbonate_orientations.push_back(
            Quaterniond(sum(0), sum(1), sum(2), sum(3)));
      }

      centroids[c] = new_centroid;
    }
  }

  return centroids;
}

void print_surface_configuration(const SurfaceState &state, int config_id) {
  trajan::log::info("=== Surface Configuration {} ===", config_id);

  for (size_t i = 0; i < state.carbonate_orientations.size(); ++i) {
    const auto &q = state.carbonate_orientations[i];

    occ::Mat3 R = q.toRotationMatrix();
    double rotation_z = std::atan2(R(1, 0), R(0, 0)) * 180.0 / M_PI;
    if (rotation_z < 0)
      rotation_z += 120.0;

    // double rotation_mod120 = std::fmod(rotation_z, 120.0);

    occ::Vec3 plane_normal = R.col(2);
    double tilt = std::acos(std::clamp(std::abs(plane_normal(2)), 0.0, 1.0)) *
                  180.0 / M_PI;

    trajan::log::info(
        "  Carbonate {}: rotation={:.1f}° (mod 120), tilt={:.1f}°", i,
        rotation_z, tilt);
  }
}

void run_qc_train(QCOpts const &opts, Trajectory &traj,
                  const Pipeline &pipeline) {
  const auto &parsed_sel = io::selection_validator(opts.carbonate_selection,
                                                   core::ATOM_RESTRICTIONS);

  trajan::log::info(
      "=== TRAINING MODE: Discovering reference configurations ===");

  // Collect trajectory states
  std::vector<SurfaceState> trajectory_states;
  int frame_num = 0;

  while (traj.next_frame()) {
    pipeline.apply(traj.frame());
    const auto &carbonates = traj.get_molecules(parsed_sel);

    SurfaceState state;
    state.frame_number = frame_num;

    for (const auto &carb : carbonates) {
      Quaterniond q = compute_carbonate_orientation(carb, traj.unit_cell());

      if (opts.use_c3_symmetry) {
        q = apply_c3_symmetry(q);
      }

      state.carbonate_orientations.push_back(q);
    }

    trajectory_states.push_back(state);
    frame_num++;
  }

  trajan::log::info("Collected {} frames with {} carbonates",
                    trajectory_states.size(),
                    trajectory_states[0].carbonate_orientations.size());

  // Cluster to find configurations
  int k = opts.n_clusters;
  trajan::log::info("Clustering into k={} configurations...", k);

  std::vector<SurfaceState> reference_configs =
      kmeans_surface_states(trajectory_states, k, opts.max_kmeans_iter);

  // Print discovered configurations
  for (int i = 0; i < k; ++i) {
    trajan::log::info("");
    print_surface_configuration(reference_configs[i], i);
  }

  // Save to file
  save_reference_configs(reference_configs, opts.reference_config_file);

  // Also show statistics from this training trajectory
  trajan::log::info("\n=== Training Trajectory Statistics ===");

  std::vector<int> assignments(trajectory_states.size());
  std::vector<int> config_counts(k, 0);

  for (size_t frame = 0; frame < trajectory_states.size(); ++frame) {
    double min_dist = 1e9;
    int best_config = 0;

    for (int c = 0; c < k; ++c) {
      double dist = surface_state_distance(trajectory_states[frame],
                                           reference_configs[c]);
      if (dist < min_dist) {
        min_dist = dist;
        best_config = c;
      }
    }

    assignments[frame] = best_config;
    config_counts[best_config]++;
  }

  for (int c = 0; c < k; ++c) {
    double percent = 100.0 * config_counts[c] / trajectory_states.size();
  }
}

void run_qc_analyse(QCOpts const &opts, Trajectory &traj,
                    const Pipeline &pipeline) {
  const auto &parsed_sel = io::selection_validator(opts.carbonate_selection,
                                                   core::ATOM_RESTRICTIONS);

  trajan::log::info(
      "=== ANALYSIS MODE: Using saved reference configurations ===");

  // LOAD reference configurations from file
  std::vector<SurfaceState> reference_configs =
      load_reference_configs(opts.reference_config_file);

  int k = reference_configs.size();

  trajan::log::info("Loaded {} reference configurations", k);
  for (int i = 0; i < k; ++i) {
    print_surface_configuration(reference_configs[i], i);
  }

  // DEBUG: Print loaded quaternions
  trajan::log::debug("Loaded quaternions from ref file:");
  for (int c = 0; c < k; ++c) {
    for (size_t i = 0; i < reference_configs[c].carbonate_orientations.size();
         ++i) {
      const auto &q = reference_configs[c].carbonate_orientations[i];
    }
  }

  trajan::log::info("\nProcessing trajectory...");

  std::vector<core::Molecule> carbonates;
  size_t frame_num = 0;

  std::vector<int> frame_assignments;
  std::vector<double> frame_distances;
  std::vector<int> config_counts(k + 1, 0);

  while (traj.next_frame()) {
    pipeline.apply(traj.frame());
    if (frame_num == 0 || traj.topology_has_changed()) {
      carbonates = traj.get_molecules(parsed_sel);
    } else {
      carbonates = traj.get_molecules(parsed_sel);
    }
    // trajan::log::info(carbonates[0].atoms()[0].repr());

    SurfaceState state;
    state.frame_number = frame_num;

    for (const auto &carb : carbonates) {
      Quaterniond q = compute_carbonate_orientation(carb, traj.unit_cell());
      if (opts.use_c3_symmetry) {
        q = apply_c3_symmetry(q);
      }
      state.carbonate_orientations.push_back(q);
    }
    // print_surface_configuration(state, 0);

    if (state.carbonate_orientations.size() !=
        reference_configs[0].carbonate_orientations.size()) {
      trajan::log::error("Frame {} has {} carbonates, but reference has {}",
                         frame_num, state.carbonate_orientations.size(),
                         reference_configs[0].carbonate_orientations.size());
      throw std::runtime_error("Carbonate count mismatch");
    }

    // print_surface_configuration(state, frame_num);

    double min_dist = 1e9;
    int best_config = 0;
    for (int c = 0; c < k; ++c) {
      double dist =
          surface_state_distance_with_permutation(state, reference_configs[c]);
      if (dist < min_dist) {
        min_dist = dist;
        best_config = c + 1;
      }
    }
    if (min_dist > opts.distance_threshold) {
      best_config = 0;
    }
    frame_assignments.push_back(best_config);
    frame_distances.push_back(min_dist);
    config_counts[best_config]++;
    frame_num++;
  }

  trajan::log::info("\n=== Configuration Occupancies ===");
  for (int c = 0; c < k + 1; ++c) {
    double percent = 100.0 * config_counts[c] / frame_assignments.size();
    if (c > 0) {
      trajan::log::info("  Config {}: {:.1f}% ({} frames)", c, percent,
                        config_counts[c]);
    } else {
      trajan::log::info("  None   {}: {:.1f}% ({} frames)", c, percent,
                        config_counts[c]);
    }
  }

  trajan::log::info("\n=== Quality Check ===");
  double avg_distance = 0.0;
  double min_distance = frame_distances[0];
  double max_distance = 0.0;
  int min_dist_frame = 0;
  int max_dist_frame = 0;

  for (size_t i = 0; i < frame_distances.size(); ++i) {
    avg_distance += frame_distances[i];
    if (frame_distances[i] > max_distance) {
      max_distance = frame_distances[i];
      max_dist_frame = i;
    }
    if (frame_distances[i] < min_distance) {
      min_distance = frame_distances[i];
      min_dist_frame = i;
    }
  }
  avg_distance /= frame_distances.size();

  trajan::log::info("Average distance to nearest config: {:.4f} rad",
                    avg_distance);
  trajan::log::info("Min distance: {:.4f} rad (frame {})", min_distance,
                    min_dist_frame);
  trajan::log::info("Max distance: {:.4f} rad (frame {})", max_distance,
                    max_dist_frame);

  if (max_distance > opts.distance_threshold) {
    trajan::log::warn(
        "Frame {} is poorly matched (distance {:.4f} > threshold {:.4f})",
        max_dist_frame, max_distance, opts.distance_threshold);
    trajan::log::warn(
        "This may indicate a novel configuration not in training set");
  }

  // Detect transitions
  trajan::log::info("\n=== Configuration Transitions ===");
  int prev_config = frame_assignments[0];
  int n_transitions = 0;

  for (size_t frame = 1; frame < frame_assignments.size(); ++frame) {
    int curr_config = frame_assignments[frame];
    if (curr_config != prev_config) {
      trajan::log::info("Frame {}: Config {} → Config {}", frame, prev_config,
                        curr_config);
      n_transitions++;
      prev_config = curr_config;
    }
  }

  trajan::log::info("Total transitions: {}", n_transitions);

  // Save results if output file specified
  if (!opts.output_file.empty()) {
    std::ofstream out(opts.output_file);
    out << "# Frame Configuration Distance\n";
    for (size_t frame = 0; frame < frame_assignments.size(); ++frame) {
      out << frame << " " << frame_assignments[frame] << " "
          << frame_distances[frame] << "\n";
    }
    out.close();
    trajan::log::info("Results saved to {}", opts.output_file);
  }
}

void run_qc_subcommand(QCOpts const &opts, Trajectory &traj,
                       const Pipeline &pipeline) {
  if (opts.train_mode) {
    run_qc_train(opts, traj, pipeline);
  } else {
    run_qc_analyse(opts, traj, pipeline);
  }
}

CLI::App *add_qc_subcommand(CLI::App &app, Trajectory &traj,
                             Pipeline &pipeline) {
  CLI::App *qc =
      app.add_subcommand("qc", "Quaternion clustering of carbonates");
  auto opts = std::make_shared<QCOpts>();
  qc->add_option("--carbonates", opts->carbonate_selection,
                 "First selection (prefix: i=atom indices, a=atom types, "
                 "j=molecule indices, m=molecule types)\n"
                 "Examples:\n"
                 "  j1,3-5      (molecule indices 1,3,4,5)\n"
                 "  mM1,M2      (molecule types M1,M2)")
      ->required();

  qc->add_flag("--train", opts->train_mode,
               "Training mode: discover and save reference configurations");

  qc->add_flag("--c3-symmetry", opts->use_c3_symmetry,
               "Apply C3 symmetry (account for 120° rotational equivalence)")
      ->default_val(true);

  qc->add_option("-k,--clusters", opts->n_clusters,
                 "Number of surface configurations to identify")
      ->default_val(3);

  qc->add_option("--ref-config", opts->reference_config_file,
                 "Reference configuration file (.json)")
      ->default_val("surface_configs.json");

  qc->add_option("-o,--output", opts->output_file,
                 "Output file for frame assignments");

  qc->add_option("--distance-threshold", opts->distance_threshold,
                 "Warn if frame distance exceeds this (radians)")
      ->default_val(0.3);

  qc->callback([opts, &traj, &pipeline]() {
    trajan::log::set_subcommand_log_pattern("qc");
    run_qc_subcommand(*opts, traj, pipeline);
  });
  return qc;
}

} // namespace trajan::main
