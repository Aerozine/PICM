if(USE_GPU AND USE_CUDA_FOUND)
    enable_language(CUDA)

    add_library(GPIC_GPU STATIC GPIC.cu)

    target_include_directories(GPIC_GPU PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/../.."
        "${CMAKE_CURRENT_SOURCE_DIR}/../../core"
        "${CMAKE_CURRENT_SOURCE_DIR}/../../core/IterativeMethod")

    target_compile_definitions(GPIC_GPU PUBLIC
        USE_CUDA
        $<IF:$<BOOL:${USE_FLOAT_PRECISION}>,USE_FLOAT,USE_DOUBLE>
        $<IF:$<BOOL:${USE_FAST_ALG}>,USE_SPEED,USE_PRECISION>)

    set_target_properties(GPIC_GPU PROPERTIES
        CUDA_ARCHITECTURES "${PICM_CUDA_ARCHITECTURES}")

    set_target_properties(GPIC_GPU PROPERTIES
        CUDA_STANDARD          17
        CUDA_STANDARD_REQUIRED ON
        CUDA_SEPARABLE_COMPILATION OFF)

    target_compile_options(GPIC_GPU PRIVATE
        $<$<CONFIG:Release>:--use_fast_math -O3>
        $<$<CONFIG:Debug>:-G -g -O0>)

    # Pull in the iterative method GPU library for MakeIncompressible dispatch.
    target_link_libraries(GPIC_GPU PUBLIC IterativeMethod_GPU)
endif()
