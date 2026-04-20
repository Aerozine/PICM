file(GLOB SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/solvers/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/solvers/*/*.cpp")

add_subdirectory(core/IterativeMethod)

add_executable(PIC ${SOURCES})

set_target_properties(PIC PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

target_include_directories(PIC PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/core"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/IterativeMethod"
    "${CMAKE_CURRENT_SOURCE_DIR}/solvers")

target_link_libraries(PIC PRIVATE
    nlohmann_json::nlohmann_json
    IterativeMethod_CPU)           # always linked

if(USE_CUDA_FOUND)
    target_link_libraries(PIC PRIVATE IterativeMethod_GPU CUDA::cudart)
    target_compile_definitions(PIC PRIVATE USE_CUDA)
endif()

if(ZLIB_FOUND)
    target_link_libraries(PIC PRIVATE ZLIB::ZLIB)
    target_compile_definitions(PIC PRIVATE HAVE_ZLIB)
endif()

if(OpenMP_CXX_FOUND)
    target_link_libraries(PIC PRIVATE OpenMP::OpenMP_CXX)
    target_compile_definitions(PIC PRIVATE USE_OPENMP)
endif()

target_compile_definitions(PIC PRIVATE
    $<IF:$<BOOL:${USE_FLOAT_PRECISION}>,USE_FLOAT,USE_DOUBLE>)

target_compile_options(PIC PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:MSVC>:/W4>)

target_compile_options(PIC PRIVATE
    $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:-O3 -march=native>
    $<$<AND:$<CONFIG:Debug>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:-O0 -g>)
