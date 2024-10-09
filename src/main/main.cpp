#include <CLI/CLI.hpp>
#include <trajan/core/log.h>
#include <trajan/main/input.h>
#include <trajan/main/trajan.h>
#include <trajan/main/version.h>

int main(int argc, char *argv[]) {
  trajan::log::setup_logging(2); // default logger

  CLI::App app(
      "trajan - A program for analysing molecular dyanmics trajectories");
  app.set_help_all_flag("--help-all", "Show help for all sub commands");

  // logging verbosity
  auto *verbosity_option = app.add_flag_function(
      "--verbosity{2}",
      [](int verbosity) { trajan::log::setup_logging(verbosity); },
      "logging verbosity {0=silent,1=minimal,2=normal,3=verbose,4=debug}");
  verbosity_option->default_val(2);
  verbosity_option->run_callback_for_default();
  verbosity_option->force_callback();

  // main object
  auto trajan = std::make_shared<trajan::TRAJAN>();

  // add all the subcommands here
  auto *inp = trajan::main::add_input_subcommand(app, trajan);

  // ensure we have a subcommand
  app.require_subcommand(/* min */ 0, /* max */ 0);

  // display default info
  trajan::main::print_header();

  constexpr auto *error_format = "exception:\n    {}\nterminating program.\n ";
  try {
    CLI11_PARSE(app, argc, argv);
  } catch (const char *ex) {
    trajan::log::error(error_format, ex);
    spdlog::dump_backtrace();
    return 1;
  } catch (std::string &ex) {
    trajan::log::error(error_format, ex);
    spdlog::dump_backtrace();
    return 1;
  } catch (std::exception &ex) {
    trajan::log::error(error_format, ex.what());
    spdlog::dump_backtrace();
    return 1;
  } catch (...) {
    trajan::log::error("Exception:\n- Unknown...\n");
    spdlog::dump_backtrace();
    return 1;
  }

  trajan::log::info("TRAJAN terminated successfully.");

  return 0;
}
