cpmaddpackage(
  NAME
  fmt
  GITHUB_REPOSITORY
  "fmtlib/fmt"
  VERSION
  "10.2.0"
  GIT_TAG
  "10.2.0")

cpmaddpackage(NAME CLI11 GITHUB_REPOSITORY "CLIUtils/CLI11" VERSION "2.4.1")

cpmaddpackage(
  NAME
  spdlog
  GITHUB_REPOSITORY
  "gabime/spdlog"
  OPTIONS
  "SPDLOG_FMT_EXTERNAL ON"
  VERSION
  "1.12.0")
