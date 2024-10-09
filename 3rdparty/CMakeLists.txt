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

cpmaddpackage(
  NAME eigen3 URL
  "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip"
  DOWNLOAD_ONLY YES # Eigen CMakelists.txt not intended for library use
)

# Add the Eigen library target to link against
if("${eigen3_ADDED}")
  add_library(Eigen3::Eigen INTERFACE IMPORTED GLOBAL)
  target_include_directories(Eigen3::Eigen INTERFACE "${eigen3_SOURCE_DIR}")
endif()
