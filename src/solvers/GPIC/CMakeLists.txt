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

    # Mirror the GPU arch handling from IterativeMethod/CMakeLists.txt.
    if(GPU_ARCH)
        string(REGEX REPLACE "^sm_" "" _arch_num "${GPU_ARCH}")
        set_target_properties(GPIC_GPU PROPERTIES
            CUDA_ARCHITECTURES "${_arch_num}")
    else()
        if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.24")
            set_target_properties(GPIC_GPU PROPERTIES
                CUDA_ARCHITECTURES "native")
        else()
            set_target_properties(GPIC_GPU PROPERTIES
                CUDA_ARCHITECTURES "70;75;80;89")
        endif()
    endif()

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
