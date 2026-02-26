cmake_minimum_required(VERSION 3.15)
project(PIC)

option(USE_FLOAT_PRECISION "Use float instead of double" OFF)

include(FetchContent)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

message(STATUS "Setting nlohmannJson")
find_package(nlohmann_json QUIET)
if(NOT TARGET nlohmann_json::nlohmann_json)
  message(STATUS "Fetching nlohmann-json")
  FetchContent_Declare(json GIT_REPOSITORY https://github.com/nlohmann/json.git GIT_TAG v3.12.0)
  FetchContent_MakeAvailable(json)
endif()
# Find OpenMP
find_package(OpenMP REQUIRED)

find_package(ZLIB QUIET)
if(ZLIB_FOUND)
  message(STATUS "ZLib found (${ZLIB_VERSION_STRING}) – VTI output will be compressed")
else()
  message(STATUS "ZLib NOT found – VTI output will use raw binary (no compression)")
endif()
add_subdirectory(src)
