add_executable(
  tacpp
  "${CMAKE_CURRENT_SOURCE_DIR}/tacpp.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/tacpp_input.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/version.cpp")
target_link_libraries(tacpp PUBLIC tacpp_core tacpp_io CLI11::CLI11)

set_target_properties(tacpp PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                       "${CMAKE_BINARY_DIR}/bin")

target_compile_features(tacpp PUBLIC cxx_std_20)
