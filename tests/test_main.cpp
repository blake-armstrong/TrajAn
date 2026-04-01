#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <fstream>
#include <occ/core/element.h>
#include <occ/crystal/unitcell.h>
#include <occ/core/units.h>
#include <trajan/core/frame.h>
#include <trajan/core/log.h>
#include <trajan/core/pipeline.h>
#include <trajan/core/trajectory.h>
#include <trajan/io/xyz.h>
#include <trajan/main/trajan_modify.h>
#include <trajan/main/trajan_write.h>

namespace fs = std::filesystem;

using Catch::Approx;
using trajan::core::Atom;
using trajan::core::Frame;
using trajan::core::FrameTransform;
using trajan::core::Pipeline;
using trajan::core::Trajectory;
using trajan::main::ModifyOpts;
using trajan::main::WriteOpts;

// ── helpers ──────────────────────────────────────────────────────────────────

static std::vector<Atom> make_water_atoms() {
  using occ::core::Element;
  std::vector<Atom> atoms;
  atoms.emplace_back(occ::Vec3{0.000, 0.000, 0.000}, Element("O", true), 0);
  atoms[0].type = "O";
  atoms.emplace_back(occ::Vec3{0.960, 0.000, 0.000}, Element("H", true), 1);
  atoms[1].type = "H";
  atoms.emplace_back(occ::Vec3{-0.240, 0.927, 0.000}, Element("H", true), 2);
  atoms[2].type = "H";
  return atoms;
}

// Write an XYZ file with n_frames frames; each frame shifts x by frame_idx Å.
static fs::path write_xyz_frames(const fs::path &path, int n_frames) {
  trajan::io::XYZHandler writer;
  writer.set_file_path(path);
  writer.initialise(trajan::io::FileHandler::Mode::Write);
  for (int f = 0; f < n_frames; ++f) {
    Frame frame;
    auto atoms = make_water_atoms();
    for (auto &atom : atoms)
      atom.x += static_cast<double>(f);
    frame.set_atoms(atoms);
    writer.write_frame(frame);
  }
  writer.finalise();
  return path;
}

// ── Pipeline unit tests ───────────────────────────────────────────────────────

TEST_CASE("Pipeline - empty pipeline is a no-op", "[unit][pipeline]") {
  Pipeline pipeline;
  REQUIRE(pipeline.empty());
  REQUIRE(pipeline.size() == 0);

  Frame frame;
  frame.set_atoms(make_water_atoms());
  double ox_before = frame.atoms()[0].x;

  pipeline.apply(frame);

  REQUIRE(frame.atoms()[0].x == Approx(ox_before));
}

TEST_CASE("Pipeline - single transform is applied", "[unit][pipeline]") {
  Pipeline pipeline;
  pipeline.add_transform([](Frame &f) {
    for (auto &atom : f.atoms())
      atom.x += 1.0;
  });
  REQUIRE(!pipeline.empty());
  REQUIRE(pipeline.size() == 1);

  Frame frame;
  frame.set_atoms(make_water_atoms());
  double ox_before = frame.atoms()[0].x;

  pipeline.apply(frame);

  REQUIRE(frame.atoms()[0].x == Approx(ox_before + 1.0));
}

TEST_CASE("Pipeline - multiple transforms are applied in order",
          "[unit][pipeline]") {
  Pipeline pipeline;

  // Transform 1: set x = 10
  pipeline.add_transform(
      [](Frame &f) { f.atoms()[0].x = 10.0; });

  // Transform 2: multiply x by 2
  pipeline.add_transform(
      [](Frame &f) { f.atoms()[0].x *= 2.0; });

  REQUIRE(pipeline.size() == 2);

  Frame frame;
  frame.set_atoms(make_water_atoms());
  pipeline.apply(frame);

  // Should be 10 * 2 = 20, not 2 * 10
  REQUIRE(frame.atoms()[0].x == Approx(20.0));
}

// ── Translate transform tests ─────────────────────────────────────────────────

TEST_CASE("register_modify_transforms - translate shifts all atoms",
          "[unit][modify]") {
  Pipeline pipeline;
  ModifyOpts opts;
  opts.translate = {1.0, 2.0, 3.0};
  opts.wrap = false;
  Trajectory traj;
  trajan::main::register_modify_transforms(opts, traj, pipeline);

  REQUIRE(pipeline.size() == 1);

  Frame frame;
  frame.set_atoms(make_water_atoms());
  const double ox0 = frame.atoms()[0].x;
  const double oy0 = frame.atoms()[0].y;
  const double oz0 = frame.atoms()[0].z;

  pipeline.apply(frame);

  REQUIRE(frame.atoms()[0].x == Approx(ox0 + 1.0));
  REQUIRE(frame.atoms()[0].y == Approx(oy0 + 2.0));
  REQUIRE(frame.atoms()[0].z == Approx(oz0 + 3.0));
  // All three atoms should be shifted
  REQUIRE(frame.atoms()[1].x == Approx(0.960 + 1.0));
  REQUIRE(frame.atoms()[2].x == Approx(-0.240 + 1.0));
}

TEST_CASE("register_modify_transforms - zero translate registers no transform",
          "[unit][modify]") {
  Pipeline pipeline;
  ModifyOpts opts;
  opts.translate = {0.0, 0.0, 0.0};
  opts.wrap = false;
  Trajectory traj;
  trajan::main::register_modify_transforms(opts, traj, pipeline);

  REQUIRE(pipeline.empty());
}

TEST_CASE("register_modify_transforms - wrap into unit cell",
          "[unit][modify]") {
  Pipeline pipeline;
  ModifyOpts opts;
  opts.translate = {0.0, 0.0, 0.0};
  opts.wrap = true;
  Trajectory traj;
  trajan::main::register_modify_transforms(opts, traj, pipeline);

  REQUIRE(pipeline.size() == 1);

  // Build a frame with a 10 Å cubic unit cell; place the O atom at 12 Å (just
  // outside the box) so wrap should bring it to 2 Å.
  Frame frame;
  auto atoms = make_water_atoms();
  atoms[0].x = 12.0; // outside [0, 10)
  frame.set_atoms(atoms);

  const double box = 10.0;
  auto uc = occ::crystal::triclinic_cell(box, box, box,
                                         occ::units::radians(90.0),
                                         occ::units::radians(90.0),
                                         occ::units::radians(90.0));
  frame.set_unit_cell(uc);

  pipeline.apply(frame);

  REQUIRE(frame.atoms()[0].x == Approx(2.0).margin(1e-4));
}

TEST_CASE("register_modify_transforms - wrap without unit cell is a no-op",
          "[unit][modify]") {
  Pipeline pipeline;
  ModifyOpts opts;
  opts.translate = {0.0, 0.0, 0.0};
  opts.wrap = true;
  Trajectory traj;
  trajan::main::register_modify_transforms(opts, traj, pipeline);

  Frame frame;
  auto atoms = make_water_atoms();
  atoms[0].x = 50.0;
  frame.set_atoms(atoms);
  // No unit cell set

  // Should not throw; atom position unchanged because wrap is skipped
  REQUIRE_NOTHROW(pipeline.apply(frame));
  REQUIRE(frame.atoms()[0].x == Approx(50.0));
}

TEST_CASE("register_modify_transforms - translate and wrap both registered",
          "[unit][modify]") {
  Pipeline pipeline;
  ModifyOpts opts;
  opts.translate = {5.0, 0.0, 0.0}; // shift into box then wrap
  opts.wrap = true;
  Trajectory traj;
  trajan::main::register_modify_transforms(opts, traj, pipeline);

  REQUIRE(pipeline.size() == 2);
}

// ── Write subcommand tests ────────────────────────────────────────────────────

TEST_CASE("run_write_subcommand - streaming: writes all frames to XYZ",
          "[unit][write]") {
  fs::path src = "/tmp/test_write_src.xyz";
  fs::path dst = "/tmp/test_write_dst.xyz";
  constexpr int n_frames = 3;

  write_xyz_frames(src, n_frames);

  {
    Trajectory traj;
    traj.load_files({src});

    Pipeline pipeline; // empty — no transforms
    WriteOpts opts;
    opts.outfile = dst;
    trajan::main::run_write_subcommand(opts, traj, pipeline);
  }

  // Read back and count frames / check positions
  trajan::io::XYZHandler reader;
  reader.set_file_path(dst);
  reader.initialise(trajan::io::FileHandler::Mode::Read);

  Frame frame;
  frame.set_atoms(make_water_atoms());

  int count = 0;
  std::vector<double> ox_values;
  while (reader.read_frame(frame)) {
    ox_values.push_back(frame.atoms()[0].x);
    count++;
  }
  reader.finalise();

  REQUIRE(count == n_frames);
  for (int f = 0; f < n_frames; ++f) {
    REQUIRE(ox_values[f] == Approx(static_cast<double>(f)).margin(1e-4));
  }
}

TEST_CASE("run_write_subcommand - in-memory: writes all frames to XYZ",
          "[unit][write]") {
  fs::path src = "/tmp/test_write_inmem_src.xyz";
  fs::path dst = "/tmp/test_write_inmem_dst.xyz";
  constexpr int n_frames = 3;

  write_xyz_frames(src, n_frames);

  {
    Trajectory traj;
    traj.load_files_into_memory({src});

    Pipeline pipeline;
    WriteOpts opts;
    opts.outfile = dst;
    trajan::main::run_write_subcommand(opts, traj, pipeline);
  }

  trajan::io::XYZHandler reader;
  reader.set_file_path(dst);
  reader.initialise(trajan::io::FileHandler::Mode::Read);

  Frame frame;
  frame.set_atoms(make_water_atoms());

  int count = 0;
  std::vector<double> ox_values;
  while (reader.read_frame(frame)) {
    ox_values.push_back(frame.atoms()[0].x);
    count++;
  }
  reader.finalise();

  REQUIRE(count == n_frames);
  for (int f = 0; f < n_frames; ++f) {
    REQUIRE(ox_values[f] == Approx(static_cast<double>(f)).margin(1e-4));
  }
}

// ── Modify + Write pipeline integration tests ─────────────────────────────────

TEST_CASE(
    "modify + write pipeline - translate is applied before writing (streaming)",
    "[unit][pipeline][modify][write]") {
  fs::path src = "/tmp/test_modify_write_src.xyz";
  fs::path dst = "/tmp/test_modify_write_dst.xyz";
  constexpr int n_frames = 2;
  constexpr double shift = 5.0;

  write_xyz_frames(src, n_frames);

  {
    Trajectory traj;
    traj.load_files({src});

    Pipeline pipeline;
    ModifyOpts modify_opts;
    modify_opts.translate = {shift, 0.0, 0.0};
    trajan::main::register_modify_transforms(modify_opts, traj, pipeline);

    WriteOpts write_opts;
    write_opts.outfile = dst;
    trajan::main::run_write_subcommand(write_opts, traj, pipeline);
  }

  trajan::io::XYZHandler reader;
  reader.set_file_path(dst);
  reader.initialise(trajan::io::FileHandler::Mode::Read);

  Frame frame;
  frame.set_atoms(make_water_atoms());

  int count = 0;
  std::vector<double> ox_values;
  while (reader.read_frame(frame)) {
    ox_values.push_back(frame.atoms()[0].x);
    count++;
  }
  reader.finalise();

  REQUIRE(count == n_frames);
  // Frame f originally has O.x == f; after translate it should be f + shift
  for (int f = 0; f < n_frames; ++f) {
    REQUIRE(ox_values[f] == Approx(static_cast<double>(f) + shift).margin(1e-4));
  }
}

TEST_CASE(
    "modify + write pipeline - translate is applied before writing (in-memory)",
    "[unit][pipeline][modify][write]") {
  fs::path src = "/tmp/test_modify_write_mem_src.xyz";
  fs::path dst = "/tmp/test_modify_write_mem_dst.xyz";
  constexpr int n_frames = 2;
  constexpr double shift = 3.0;

  write_xyz_frames(src, n_frames);

  {
    Trajectory traj;
    traj.load_files_into_memory({src});

    Pipeline pipeline;
    ModifyOpts modify_opts;
    modify_opts.translate = {shift, 0.0, 0.0};
    trajan::main::register_modify_transforms(modify_opts, traj, pipeline);

    WriteOpts write_opts;
    write_opts.outfile = dst;
    trajan::main::run_write_subcommand(write_opts, traj, pipeline);
  }

  trajan::io::XYZHandler reader;
  reader.set_file_path(dst);
  reader.initialise(trajan::io::FileHandler::Mode::Read);

  Frame frame;
  frame.set_atoms(make_water_atoms());

  int count = 0;
  std::vector<double> ox_values;
  while (reader.read_frame(frame)) {
    ox_values.push_back(frame.atoms()[0].x);
    count++;
  }
  reader.finalise();

  REQUIRE(count == n_frames);
  for (int f = 0; f < n_frames; ++f) {
    REQUIRE(ox_values[f] == Approx(static_cast<double>(f) + shift).margin(1e-4));
  }
}

TEST_CASE("modify + write pipeline - streaming and in-memory produce identical "
          "output",
          "[unit][pipeline][modify][write]") {
  fs::path src = "/tmp/test_modify_write_parity_src.xyz";
  fs::path dst_stream = "/tmp/test_modify_write_parity_stream.xyz";
  fs::path dst_mem = "/tmp/test_modify_write_parity_mem.xyz";
  constexpr int n_frames = 4;

  write_xyz_frames(src, n_frames);

  auto run = [&](bool into_memory, const fs::path &dst) {
    Trajectory traj;
    if (into_memory)
      traj.load_files_into_memory({src});
    else
      traj.load_files({src});

    Pipeline pipeline;
    ModifyOpts modify_opts;
    modify_opts.translate = {1.5, -0.5, 2.0};
    trajan::main::register_modify_transforms(modify_opts, traj, pipeline);

    WriteOpts write_opts;
    write_opts.outfile = dst;
    trajan::main::run_write_subcommand(write_opts, traj, pipeline);
  };

  run(false, dst_stream);
  run(true, dst_mem);

  // Compare the two output files line-by-line
  std::ifstream f1(dst_stream), f2(dst_mem);
  REQUIRE(f1.is_open());
  REQUIRE(f2.is_open());

  std::string line1, line2;
  int line_num = 0;
  while (std::getline(f1, line1) && std::getline(f2, line2)) {
    INFO("Line " << line_num << ": '" << line1 << "' vs '" << line2 << "'");
    REQUIRE(line1 == line2);
    line_num++;
  }
  REQUIRE(line_num > 0);
  // Both should be fully consumed at this point
  bool f1_done = !std::getline(f1, line1);
  bool f2_done = !std::getline(f2, line2);
  REQUIRE(f1_done);
  REQUIRE(f2_done);
}
