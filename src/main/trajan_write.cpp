#include <trajan/core/log.h>
#include <trajan/main/trajan_write.h>

namespace trajan::main {

void run_write_subcommand(const WriteOpts &opts, Trajectory &traj,
                          const Pipeline &pipeline) {
  traj.set_output_file(opts.outfile);
  traj.set_output_original_ids(opts.original_ids);

  size_t frame_count = 0;
  while (traj.next_frame()) {
    pipeline.apply(traj.frame(), "write");
    traj.write_frame();
    frame_count++;
  }

  trajan::log::info("Wrote {} frame(s) to {}", frame_count,
                    opts.outfile.string());
}

CLI::App *add_write_subcommand(CLI::App &app, Trajectory &traj,
                               Pipeline &pipeline) {
  CLI::App *write = app.add_subcommand(
      "write", "Write trajectory frames (with any queued modifications) to an "
               "output file. Requires load subcommand.");

  auto opts = std::make_shared<WriteOpts>();
  write->add_option("outfile", opts->outfile, "Output trajectory file path")
      ->required();
  write->add_flag("--original-ids", opts->original_ids,
                  "Use the original atom serial numbers from the input file "
                  "instead of renumbering 1..N.");

  write->callback([opts, &traj, &pipeline]() {
    trajan::log::set_subcommand_log_pattern("write");
    trajan::log::debug("Beginning write subcommand");
    run_write_subcommand(*opts, traj, pipeline);
    trajan::log::debug("write subcommand completed successfully");
  });

  return write;
}

} // namespace trajan::main
