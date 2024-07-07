add_library(
  tacpp_core "${CMAKE_CURRENT_SOURCE_DIR}/radial_distribution_func.cpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/log.cpp")

target_link_libraries(tacpp_core PUBLIC fmt::fmt spdlog::spdlog Eigen3::Eigen)

target_include_directories(tacpp_core PUBLIC "${TACPP_INCLUDE_DIR}")

target_compile_features(tacpp_core PUBLIC cxx_std_20)
