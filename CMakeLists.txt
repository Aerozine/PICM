cmake_minimum_required(VERSION 3.15)
project(PIC)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

message(STATUS "Setting up Eigen")
include(FetchContent)
find_package(Eigen3 5.0 QUIET)
if(NOT TARGET Eigen3::Eigen)
  if(EXISTS "/usr/include/eigen3/Eigen/Core")
    message(STATUS "Using system Eigen headers")
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(Eigen3::Eigen PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                                   /usr/include/eigen3)
  else()
    message(STATUS "Fetching Eigen")
    FetchContent_Declare(
      eigen
      GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
      GIT_TAG 5.0.1)
    FetchContent_MakeAvailable(eigen)
  endif()
endif()
message(STATUS "Setting nlohmannJson")
include(FetchContent)

find_package(nlohmann_json QUIET)
if(NOT TARGET nlohmann_json::nlohmann_json)
  message(STATUS "Fetching nlohmann-json")
  FetchContent_Declare(
    json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG v3.12.0 # latest version
  )
  FetchContent_MakeAvailable(json)
endif()
# KOKKOS ? find_package(Kokkos CONFIG) # Try to find Kokkos externally
# if(Kokkos_FOUND) message(STATUS "Found Kokkos: ${Kokkos_DIR} (version
# \"${Kokkos_VERSION}\")") else() message(STATUS "Kokkos not found externally.
# Fetching via FetchContent.") include(FetchContent) FetchContent_Declare(
# Kokkos URL https://github.com/kokkos/kokkos/archive/refs/tags/4.4.01.tar.gz )
# FetchContent_MakeAvailable(Kokkos) endif() ...
add_subdirectory(src)
