cmake_minimum_required(VERSION 3.16.0 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(PROJECT_VERSION_MAJOR "0")
set(PROJECT_VERSION_MINOR "0")
set(PROJECT_VERSION_PATCH "0")

project(
  trajan
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
  option(WITH_PYTHON_BINDINGS "Build python bindings" OFF)
endif()

if("${WITH_PYTHON_BINDINGS}")
    find_package(Python REQUIRED COMPONENTS Interpreter Development.Module
	OPTIONAL_COMPONENTS Development.SABIModule)
    message(STATUS "Enabled python bindings, setting CMAKE_POSITION_INDEPENDENT_CODE globally")
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/CPM.cmake")

# TBB's GNU.cmake guards -fcf-protection=full with `aarch64` but macOS ARM
# reports `arm64` to CMake.  With GCC this flag is a hard error, so remap the
# processor name before any dependency pulls in TBB.
if(APPLE AND CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64" AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(STATUS "Remapping CMAKE_SYSTEM_PROCESSOR arm64 -> aarch64 for GCC/TBB compatibility")
    set(CMAKE_SYSTEM_PROCESSOR "aarch64")
endif()

# don't use install commands from third party
add_subdirectory(3rdparty EXCLUDE_FROM_ALL)
add_subdirectory(include)
add_subdirectory(src)

# option(BUILD_TESTING "Build the testing tree" ON) 
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
  enable_testing()
  add_subdirectory(tests)
endif()

if("${WITH_CLANG_TIDY}")
  set(CMAKE_CXX_CLANG_TIDY
      "clang-tidy;clang-diagnostic-*,clang-analyzer-*,-*,bugprone*,modernize*,performance*,-modernize-pass-by-value,-modernize-use-auto,-modernize-use-using"
  )
endif()
