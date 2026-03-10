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
  FetchContent_Declare(
    json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG v3.12.0)
  FetchContent_MakeAvailable(json)
endif()

find_package(OpenMP QUIET)
if(NOT OpenMP_CXX_FOUND AND APPLE)
  # On macOS , libomp is available via Homebrew (both intel and Apple-silicon
  # paths are checked automatically)
  find_program(BREW brew)
  if(BREW)
    execute_process(
      COMMAND ${BREW} --prefix libomp
      OUTPUT_VARIABLE LIBOMP_PREFIX
      OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  endif()

  if(EXISTS "${LIBOMP_PREFIX}/include/omp.h")
    message(STATUS "OpenMP: using Homebrew libomp at ${LIBOMP_PREFIX}")
    add_library(OpenMP::OpenMP_CXX INTERFACE IMPORTED)
    target_compile_options(
      OpenMP::OpenMP_CXX INTERFACE -Xpreprocessor -fopenmp
                                   -I${LIBOMP_PREFIX}/include)
    target_link_libraries(OpenMP::OpenMP_CXX INTERFACE -L${LIBOMP_PREFIX}/lib
                                                       -lomp)
    set(OpenMP_CXX_FOUND TRUE)
  endif()
endif()

if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP: enabled")
else()
  message(STATUS "OpenMP: NOT found - building single-threaded")
endif()

find_package(ZLIB QUIET)
if(ZLIB_FOUND)
  message(
    STATUS "ZLib found (${ZLIB_VERSION_STRING}) – VTI output will be compressed"
  )
else()
  message(
    STATUS "ZLib NOT found – VTI output will use raw binary (no compression)")
endif()
add_subdirectory(src)
