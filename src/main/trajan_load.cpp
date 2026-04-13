#include <trajan/main/trajan_load.h>
#include <trajan/core/log.h>
#include <stdexcept>
#include <string>

namespace trajan::main {

// Parse a Python-style slice string "b:e:s" into (begin, end, step).
// Missing components use defaults: begin=0, end=-1 (∞), step=1.
static void parse_slice(const std::string &s, int &begin, int &end, int &step) {
  begin = 0;
  end = -1;
  step = 1;
  if (s.empty())
    return;

  // Split on ':'
  std::vector<std::string> parts;
  std::string cur;
  for (char c : s) {
    if (c == ':') {
      parts.push_back(cur);
      cur.clear();
    } else {
      cur += c;
    }
  }
  parts.push_back(cur);

  if (parts.size() > 3)
    throw std::invalid_argument("--slice: too many ':' separators (expected b:e:s)");

  auto parse_int = [](const std::string &tok, int def) -> int {
    if (tok.empty()) return def;
    return std::stoi(tok);
  };

  begin = parse_int(parts[0], 0);
  if (parts.size() >= 2) end  = parse_int(parts[1], -1);
  if (parts.size() >= 3) step = parse_int(parts[2], 1);

  if (step <= 0)
    throw std::invalid_argument("--slice: step must be >= 1");
  if (begin < 0)
    throw std::invalid_argument("--slice: begin must be >= 0");
}

void run_load_subcommand(const LoadOpts &opts, Trajectory &traj) {
  if (opts.into_mem) {
    trajan::log::debug("Files will be loaded into memory");
    traj.load_files_into_memory(opts.infiles);
  } else {
    trajan::log::debug("Files will NOT be loaded into memory.");
    trajan::log::debug("File paths saved for analysis in next subcommand.");
    traj.load_files(opts.infiles);
  }

  if (!opts.slice.empty()) {
    int b, e, s;
    parse_slice(opts.slice, b, e, s);
    trajan::log::info("load: frame slice [{} : {} : {}] (end -1 = no limit)",
                      b, e, s);
    traj.set_slice(b, e, s);
  }
}

CLI::App *add_load_subcommand(CLI::App &app, Trajectory &traj) {
  CLI::App *load =
      app.add_subcommand("load", "Load trajectory data into program for "
                                 "analysis. Required by most subcommands.");
  auto opts = std::make_shared<LoadOpts>();
  load->add_option("files", opts->infiles, "Input trajectory file names")
      ->check(CLI::ExistingFile);
  load->add_flag("--into-mem", opts->into_mem,
                 "Whether or not to load files into memory. This is not "
                 "recommended for command-line use.");
  load->add_option("--slice", opts->slice,
                   "Frame slice in Python notation: b:e:s (begin, end, step). "
                   "Frames are 0-indexed; end is exclusive; any component may "
                   "be omitted. Examples: '1:' skips frame 0, '::2' takes "
                   "every other frame, '1:10:2' takes odd frames 1-9.");
  load->callback([opts, &traj]() {
    trajan::log::set_subcommand_log_pattern("load");
    trajan::log::debug("Beginning load subcommand");
    run_load_subcommand(*opts, traj);
    trajan::log::debug("load subcommand completed successfully");
  });
  return load;
}

} // namespace trajan::main
