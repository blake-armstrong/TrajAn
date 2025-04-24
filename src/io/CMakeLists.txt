add_library(trajan_io "${CMAKE_CURRENT_SOURCE_DIR}/pdb.cpp"
                     "${CMAKE_CURRENT_SOURCE_DIR}/dcd.cpp"
                     "${CMAKE_CURRENT_SOURCE_DIR}/file.cpp"
                     "${CMAKE_CURRENT_SOURCE_DIR}/selection.cpp"
)

target_link_libraries(trajan_io PUBLIC Eigen3::Eigen fmt::fmt spdlog::spdlog OpenMP::OpenMP_CXX unordered_dense::unordered_dense)

target_include_directories(trajan_io PUBLIC "${TRAJAN_INCLUDE_DIR}")

target_compile_features(trajan_io PUBLIC cxx_std_20)
