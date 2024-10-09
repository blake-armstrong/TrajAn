add_library(
  trajan_core
  "${CMAKE_CURRENT_SOURCE_DIR}/radial_distribution_func.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/log.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/atom.cpp")

target_link_libraries(trajan_core PUBLIC fmt::fmt spdlog::spdlog Eigen3::Eigen)

target_include_directories(trajan_core PUBLIC "${TRAJAN_INCLUDE_DIR}")

target_compile_features(trajan_core PUBLIC cxx_std_20)
