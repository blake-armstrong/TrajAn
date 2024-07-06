add_library(tacpp_io "${CMAKE_CURRENT_SOURCE_DIR}/read_pdb.cpp")
target_include_directories(tacpp_io PUBLIC "${TACPP_INCLUDE_DIR}")

target_compile_features(tacpp_io PUBLIC cxx_std_20)
