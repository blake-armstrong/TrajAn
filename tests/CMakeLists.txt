cpmaddpackage(NAME catch GITHUB_REPOSITORY "catchorg/Catch2" GIT_TAG v3.5.2)

list(APPEND CMAKE_MODULE_PATH "${catch_SOURCE_DIR}/extras")
include(CTest)
include(Catch)
