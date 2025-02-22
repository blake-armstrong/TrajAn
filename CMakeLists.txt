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

# option(USE_OPENMP "Use OpenMP for threading" ON)
find_package(OpenMP)


include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/CPM.cmake")
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
