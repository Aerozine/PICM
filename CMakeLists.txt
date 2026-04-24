cmake_minimum_required(VERSION 3.18)
project(PIC)

option(USE_FLOAT_PRECISION "Use float instead of double"       ON)
option(USE_FAST_ALG "prioritize speed over precision"          OFF)
option(USE_PARALLEL        "Enable OpenMP CPU parallelism"     ON)
option(USE_GPU             "Enable CUDA GPU pressure solvers"  OFF)

# Optional: user can set -DGPU_ARCH=89 or -DGPU_ARCH=sm_89 on the command line.
# Leave empty for CMake native auto-detect.
set(GPU_ARCH "" CACHE STRING
    "CUDA architecture number or sm_XX string (e.g. 89 or sm_89). \
Empty = native auto-detect.")

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# In Debug, default heavy options to OFF unless explicitly set.
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    if(NOT DEFINED USE_PARALLEL)
        set(USE_PARALLEL OFF CACHE BOOL "Enable OpenMP CPU parallelism" FORCE)
    endif()
    if(NOT DEFINED USE_GPU)
        set(USE_GPU OFF CACHE BOOL "Enable CUDA GPU pressure solvers" FORCE)
    endif()
endif()

set(CMAKE_NINJA_FORCE_RESPONSE_FILE "ON" CACHE BOOL "" FORCE)

find_package(nlohmann_json QUIET)
if(NOT TARGET nlohmann_json::nlohmann_json)
    message(STATUS "nlohmann/json: not found locally – fetching v3.12.0")
    include(FetchContent)
    FetchContent_Declare(json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG        v3.12.0)
    FetchContent_MakeAvailable(json)
endif()

set(OpenMP_CXX_FOUND FALSE)
if(USE_PARALLEL)
    find_package(OpenMP QUIET)
    # macOS Homebrew libomp fallback.
    if(NOT OpenMP_CXX_FOUND AND APPLE)
        find_program(BREW brew)
        if(BREW)
            execute_process(COMMAND ${BREW} --prefix libomp
                OUTPUT_VARIABLE LIBOMP_PREFIX
                OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
        endif()
        if(EXISTS "${LIBOMP_PREFIX}/include/omp.h")
            message(STATUS "OpenMP: using Homebrew libomp at ${LIBOMP_PREFIX}")
            add_library(OpenMP::OpenMP_CXX INTERFACE IMPORTED)
            target_compile_options(OpenMP::OpenMP_CXX INTERFACE
                -Xpreprocessor -fopenmp -I${LIBOMP_PREFIX}/include)
            target_link_libraries(OpenMP::OpenMP_CXX INTERFACE
                -L${LIBOMP_PREFIX}/lib -lomp)
            set(OpenMP_CXX_FOUND TRUE)
        endif()
    endif()

    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP: enabled")
    else()
        message(WARNING "OpenMP: NOT found – building single-threaded")
    endif()
else()
    message(STATUS "OpenMP: disabled (USE_PARALLEL=OFF)")
endif()

find_package(ZLIB QUIET)
if(ZLIB_FOUND)
    message(STATUS "ZLib: found ${ZLIB_VERSION_STRING} – VTI output compressed")
else()
    message(STATUS "ZLib: not found – VTI output will use raw binary")
endif()

set(USE_CUDA_FOUND FALSE)
if(USE_GPU)
    # Let CMake find the toolkit itself. The active WSL compiler can be
    # selected with PATH, CUDACXX, CMAKE_CUDA_COMPILER, or CUDAToolkit_ROOT.
    # Avoid hardcoding /usr/local/cuda so alternate installations (for example
    # the NVIDIA HPC SDK) can be selected without editing the project.
    find_package(CUDAToolkit QUIET)
    if(CUDAToolkit_FOUND)
        set(USE_CUDA_FOUND TRUE)
        message(STATUS "CUDA: found ${CUDAToolkit_VERSION} – GPU solvers enabled")
        message(STATUS "CUDA toolkit root: ${CUDAToolkit_ROOT}")
    else()
        message(WARNING
            "CUDA: CUDAToolkit not found – falling back to CPU solvers.\n"
            "  Hint: set CUDAToolkit_ROOT to your CUDA install, e.g.:\n"
            "    cmake -DCUDAToolkit_ROOT=/usr/local/cuda ..")
    endif()
else()
    message(STATUS "CUDA: disabled (USE_GPU=OFF)")
endif()

add_subdirectory(src)
