cpmaddpackage(NAME catch GITHUB_REPOSITORY "catchorg/Catch2" GIT_TAG v3.5.2)

list(APPEND CMAKE_MODULE_PATH "${catch_SOURCE_DIR}/extras")
include(CTest)
include(Catch)

add_executable(test_io "${CMAKE_CURRENT_SOURCE_DIR}/test_io.cpp")
target_link_libraries(test_io PRIVATE trajan_core trajan_io Catch2::Catch2WithMain)
catch_discover_tests(test_io)

add_executable(test_core "${CMAKE_CURRENT_SOURCE_DIR}/test_core.cpp" )
target_link_libraries(test_core PRIVATE trajan_core trajan_io Catch2::Catch2WithMain OpenMP::OpenMP_CXX)
catch_discover_tests(test_core)
