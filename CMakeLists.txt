cmake_minimum_required(VERSION 3.16.0 FATAL_ERROR)

set(PROJECT_VERSION_MAJOR "0")
set(PROJECT_VERSION_MINOR "0")
set(PROJECT_VERSION_PATCH "0")

project(
  tacpp
  VERSION
    "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}"
  LANGUAGES CXX)

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds are not supported. Please make a new directory (i.e. build directory) and run cmake from there."
  )
endif()

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  include(CTest)
  option(BUILD_FPIC "Build with position independent code" OFF)
  option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
  option(WITH_PYTHON_BINDINGS "Build python bindings" OFF)
  option(WITH_BENCHMARKS "Enable benchmarking tests" OFF)
  option(WITH_CLANG_TIDY "Enable clang tidy" OFF)
endif()

include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/CPM.cmake")
# don't use install commands from third party
add_subdirectory(3rdparty EXCLUDE_FROM_ALL)
add_subdirectory(include)
add_subdirectory(src)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
  enable_testing()
  add_subdirectory(tests)
endif()

if("${WITH_CLANG_TIDY}")
  set(CMAKE_CXX_CLANG_TIDY
      "clang-tidy;clang-diagnostic-*,clang-analyzer-*,-*,bugprone*,modernize*,performance*,-modernize-pass-by-value,-modernize-use-auto,-modernize-use-using"
  )
endif()
