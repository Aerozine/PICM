add_library(IterativeMethod_CPU STATIC
        RBGS.cpp
        CG.cpp
        Jacobi-GS.cpp
)

target_include_directories(IterativeMethod_CPU PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/..")

if(OpenMP_CXX_FOUND)
    target_link_libraries(IterativeMethod_CPU PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(IterativeMethod_CPU PUBLIC USE_OPENMP)
endif()

target_compile_definitions(IterativeMethod_CPU PUBLIC
        $<IF:$<BOOL:${USE_FLOAT_PRECISION}>,USE_FLOAT,USE_DOUBLE>)

target_compile_options(IterativeMethod_CPU PRIVATE
        $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:-O3 -march=native>
        $<$<AND:$<CONFIG:Debug>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:-O0 -g>
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>
        $<$<CXX_COMPILER_ID:MSVC>:/W4>)

# ── GPU library (CUDA) ────────────────────────────────────────────────────────
if(USE_CUDA_FOUND)
    enable_language(CUDA)

    add_library(IterativeMethod_GPU STATIC
            RBGS.cu
            CG.cu)

    target_include_directories(IterativeMethod_GPU PUBLIC
            "${CMAKE_CURRENT_SOURCE_DIR}"
            "${CMAKE_CURRENT_SOURCE_DIR}/..")

    target_compile_definitions(IterativeMethod_GPU PUBLIC
            USE_CUDA
            $<IF:$<BOOL:${USE_FLOAT_PRECISION}>,USE_FLOAT,USE_DOUBLE>)

    # Architecture: user-supplied string (e.g. "sm_89"), numeric shorthand
    # (e.g. "89"), or empty → CMake native auto-detect.
    if(GPU_ARCH)
        # Strip the "sm_" prefix if the user supplied it so both "89" and
        # "sm_89" are accepted.
        string(REGEX REPLACE "^sm_" "" _arch_num "${GPU_ARCH}")
        set_target_properties(IterativeMethod_GPU PROPERTIES
                CUDA_ARCHITECTURES "${_arch_num}")
    else()
        # "native" tells CMake to query the installed GPU at configure time.
        # Requires CMake >= 3.24; fall back to a safe default for older CMake.
        if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.24")
            set_target_properties(IterativeMethod_GPU PROPERTIES
                    CUDA_ARCHITECTURES "native")
        else()
            # Broadest useful default: Volta + Turing + Ampere + Ada
            set_target_properties(IterativeMethod_GPU PROPERTIES
                    CUDA_ARCHITECTURES "70;75;80;89")
        endif()
    endif()

    set_target_properties(IterativeMethod_GPU PROPERTIES
            CUDA_STANDARD          17
            CUDA_STANDARD_REQUIRED ON
            CUDA_SEPARABLE_COMPILATION OFF)

    target_compile_options(IterativeMethod_GPU PRIVATE
            $<$<CONFIG:Release>:--use_fast_math -O3>
            $<$<CONFIG:Debug>:-G -g -O0>)
endif()
