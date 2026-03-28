#include <catch2/catch_approx.hpp>
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <atomic>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/wait.h>
#include <vector>

#ifndef TRAJAN_BINARY
#error "TRAJAN_BINARY must be defined via CMake target_compile_definitions"
#endif

namespace fs = std::filesystem;

static std::string g_test_data_path;

std::string get_test_data_path() { return g_test_data_path; }

bool has_test_data() {
  return !g_test_data_path.empty() && fs::exists(g_test_data_path);
}

int main(int argc, char *argv[]) {
  Catch::Session session;

  auto cli = session.cli() |
             Catch::Clara::Opt(g_test_data_path, "path")["-d"]["--data-path"](
                 "Path to trajectory test data directory");

  session.cli(cli);

  int result = session.applyCommandLine(argc, argv);
  if (result != 0)
    return result;

  if (!g_test_data_path.empty() && !fs::exists(g_test_data_path)) {
    g_test_data_path.clear();
  }

  return session.run();
}

// ── helpers ───────────────────────────────────────────────────────────────────

struct RunResult {
  int exit_code;
  std::string stdout_output;
  std::string stderr_output;
};

static RunResult run_cmd(const std::string &cmd) {
  // Use a unique suffix per call to avoid collisions between parallel tests.
  static std::atomic<int> counter{0};
  int id = counter.fetch_add(1);
  fs::path tmp_out = fs::temp_directory_path() /
                     ("trajan_test_stdout_" + std::to_string(id) + ".txt");
  fs::path tmp_err = fs::temp_directory_path() /
                     ("trajan_test_stderr_" + std::to_string(id) + ".txt");

  std::string full_cmd =
      cmd + " >" + tmp_out.string() + " 2>" + tmp_err.string();
  int ret = std::system(full_cmd.c_str()); // NOLINT(cert-env33-c)

  auto read_file = [](const fs::path &p) -> std::string {
    std::ifstream f(p);
    if (!f.is_open())
      return {};
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
  };

  RunResult r;
  r.exit_code = WIFEXITED(ret) ? WEXITSTATUS(ret) : -1;
  r.stdout_output = read_file(tmp_out);
  r.stderr_output = read_file(tmp_err);

  fs::remove(tmp_out);
  fs::remove(tmp_err);
  return r;
}

static std::string trajan_bin() { return std::string(TRAJAN_BINARY); }

// Build a command string with shell-quoted arguments.
static std::string build_cmd(const std::vector<std::string> &args) {
  std::string cmd = "\"" + trajan_bin() + "\"";
  for (const auto &a : args) {
    cmd += " \"" + a + "\"";
  }
  return cmd;
}

// Count ATOM/HETATM lines in a PDB file.
static int count_pdb_atoms(const fs::path &p) {
  std::ifstream f(p);
  int n = 0;
  std::string line;
  while (std::getline(f, line)) {
    if (line.size() >= 4 &&
        (line.substr(0, 4) == "ATOM" || line.substr(0, 6) == "HETATM"))
      ++n;
  }
  return n;
}

// ── tests that don't need test data ──────────────────────────────────────────

TEST_CASE("CLI - no subcommand exits non-zero", "[cli][unit]") {
  auto r = run_cmd(build_cmd({}));
  CHECK(r.exit_code != 0);
}

TEST_CASE("CLI - --help exits zero", "[cli][unit]") {
  auto r = run_cmd(build_cmd({"--help"}));
  CHECK(r.exit_code == 0);
}

TEST_CASE("CLI - write without load rejected", "[cli][unit]") {
  // write requires load; CLI11 should refuse and exit non-zero
  auto r = run_cmd(build_cmd({"write", "/tmp/trajan_should_not_exist.pdb"}));
  CHECK(r.exit_code != 0);
}

// ── tests that require test data ─────────────────────────────────────────────

TEST_CASE("CLI - load + write round-trip (PDB)", "[cli][integration][pdb]") {
  if (!has_test_data()) {
    SKIP("No test data. Use --data-path <path>");
  }

  fs::path input = fs::path(get_test_data_path()) / "i00000000.pdb";
  if (!fs::exists(input)) {
    SKIP("i00000000.pdb not found in test data directory");
  }

  fs::path output = fs::temp_directory_path() / "trajan_cli_roundtrip_0.pdb";
  fs::remove(output); // ensure clean state

  auto r = run_cmd(build_cmd({"load", input.string(),
                               "write", output.string()}));

  REQUIRE(r.exit_code == 0);
  REQUIRE(fs::exists(output));

  int n_in = count_pdb_atoms(input);
  int n_out = count_pdb_atoms(output);
  CHECK(n_in > 0);
  CHECK(n_out == n_in);

  fs::remove(output);
}

TEST_CASE("CLI - load + write round-trip with --into-mem (PDB)",
          "[cli][integration][pdb]") {
  if (!has_test_data()) {
    SKIP("No test data. Use --data-path <path>");
  }

  fs::path input = fs::path(get_test_data_path()) / "i00000000.pdb";
  if (!fs::exists(input)) {
    SKIP("i00000000.pdb not found in test data directory");
  }

  fs::path out_stream = fs::temp_directory_path() / "trajan_cli_stream_1.pdb";
  fs::path out_mem = fs::temp_directory_path() / "trajan_cli_mem_1.pdb";
  fs::remove(out_stream);
  fs::remove(out_mem);

  // Streaming
  auto r1 = run_cmd(build_cmd({"load", input.string(),
                                "write", out_stream.string()}));
  REQUIRE(r1.exit_code == 0);

  // In-memory
  auto r2 = run_cmd(build_cmd({"load", "--into-mem", input.string(),
                                "write", out_mem.string()}));
  REQUIRE(r2.exit_code == 0);

  CHECK(count_pdb_atoms(out_stream) == count_pdb_atoms(out_mem));

  fs::remove(out_stream);
  fs::remove(out_mem);
}

TEST_CASE("CLI - load + modify --translate + write shifts atoms",
          "[cli][integration][pdb]") {
  if (!has_test_data()) {
    SKIP("No test data. Use --data-path <path>");
  }

  fs::path input = fs::path(get_test_data_path()) / "i00000000.pdb";
  if (!fs::exists(input)) {
    SKIP("i00000000.pdb not found in test data directory");
  }

  fs::path out_orig = fs::temp_directory_path() / "trajan_cli_orig_2.pdb";
  fs::path out_shifted = fs::temp_directory_path() / "trajan_cli_shifted_2.pdb";
  fs::remove(out_orig);
  fs::remove(out_shifted);

  // Write original
  auto r1 = run_cmd(build_cmd({"load", input.string(),
                                "write", out_orig.string()}));
  REQUIRE(r1.exit_code == 0);

  // Write shifted by (10, 0, 0)
  auto r2 = run_cmd(build_cmd({"load", input.string(),
                                "modify", "--translate", "10.0", "0.0", "0.0",
                                "write", out_shifted.string()}));
  REQUIRE(r2.exit_code == 0);

  // Parse first ATOM x-coordinate from both files and compare.
  auto get_first_x = [](const fs::path &p) -> double {
    std::ifstream f(p);
    std::string line;
    while (std::getline(f, line)) {
      if (line.size() >= 4 &&
          (line.substr(0, 4) == "ATOM" || line.substr(0, 6) == "HETATM")) {
        // PDB x coordinate is at columns 31-38 (1-based), i.e. substr(30,8).
        try {
          return std::stod(line.substr(30, 8));
        } catch (...) {
          return 0.0;
        }
      }
    }
    return 0.0;
  };

  double x_orig = get_first_x(out_orig);
  double x_shifted = get_first_x(out_shifted);
  CHECK(x_shifted - x_orig == Catch::Approx(10.0).epsilon(1e-3));

  fs::remove(out_orig);
  fs::remove(out_shifted);
}

TEST_CASE("CLI - modify without load rejected", "[cli][unit]") {
  auto r = run_cmd(build_cmd({"modify", "--translate", "1.0", "0.0", "0.0",
                               "write", "/tmp/trajan_no_load.pdb"}));
  CHECK(r.exit_code != 0);
}
