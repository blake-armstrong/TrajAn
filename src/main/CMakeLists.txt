add_executable(
  trajan
  "${CMAKE_CURRENT_SOURCE_DIR}/trajan.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/trajan_rdf.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/version.cpp")
target_link_libraries(trajan PUBLIC trajan_core trajan_io CLI11::CLI11)

set_target_properties(trajan PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                       "${CMAKE_BINARY_DIR}/bin")

target_compile_features(trajan PUBLIC cxx_std_20)
