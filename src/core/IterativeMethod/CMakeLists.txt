add_library(IterativeMethod_CPU STATIC
        RBGS.cpp
        CG.cpp
        Jacobi-GS.cpp
)
target_link_libraries(IterativeMethod_CPU PUBLIC nlohmann_json::nlohmann_json)

target_include_directories(IterativeMethod_CPU PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/..")

if(OpenMP_CXX_FOUND)
    target_link_libraries(IterativeMethod_CPU PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(IterativeMethod_CPU PUBLIC USE_OPENMP)
endif()

target_compile_definitions(IterativeMethod_CPU PUBLIC
        $<IF:$<BOOL:${USE_FLOAT_PRECISION}>,USE_FLOAT,USE_DOUBLE>)

find_package(Eigen3 QUIET NO_MODULE)
if(TARGET Eigen3::Eigen)
    target_link_libraries(IterativeMethod_CPU PUBLIC Eigen3::Eigen)
    target_compile_definitions(IterativeMethod_CPU PUBLIC HAVE_EIGEN)
    message(STATUS "Eigen3: enabled for CPU CG/MICCG0")
else()
    message(STATUS "Eigen3: not found - CPU CG/MICCG0 will fall back to RBGS")
endif()

target_compile_options(IterativeMethod_CPU PRIVATE
        $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:-O3 -march=native>
        $<$<AND:$<CONFIG:Debug>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:-O0 -g>
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>
        $<$<CXX_COMPILER_ID:MSVC>:/W4>)

if(USE_GPU)
# GPU library (CUDA)
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

    set_target_properties(IterativeMethod_GPU PROPERTIES
            CUDA_ARCHITECTURES "${PICM_CUDA_ARCHITECTURES}")

    set_target_properties(IterativeMethod_GPU PROPERTIES
            CUDA_STANDARD          17
            CUDA_STANDARD_REQUIRED ON
            CUDA_SEPARABLE_COMPILATION OFF)

    target_compile_options(IterativeMethod_GPU PRIVATE
            $<$<CONFIG:Release>:--use_fast_math -O3>
            $<$<CONFIG:Debug>:-G -g -O0>)
endif()
endif()
