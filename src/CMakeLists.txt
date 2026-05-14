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

#option(USE_LIKWID "Enable LIKWID performance markers" OFF)
#if(USE_LIKWID)
    #find_library(LIKWID_LIB likwid HINTS /usr/local/lib /usr/lib)
    #find_path(LIKWID_INCLUDE likwid-marker.h
            #HINTS /usr/local/include /usr/include)
    #if(LIKWID_LIB AND LIKWID_INCLUDE)
        #message(STATUS "LIKWID: found at ${LIKWID_LIB}")
        #target_compile_definitions(PIC PRIVATE LIKWID_PERFMON)
        #target_include_directories(PIC PRIVATE ${LIKWID_INCLUDE})
        #target_link_libraries(PIC PRIVATE ${LIKWID_LIB})
    #else()
        #message(WARNING "LIKWID: not found - markers disabled")
    #endif()
#endif()

target_link_libraries(PIC PUBLIC 
    nlohmann_json::nlohmann_json
    IterativeMethod_CPU)           # always linked
if(USE_GPU)
if(USE_CUDA_FOUND)
    target_link_libraries(PIC PRIVATE IterativeMethod_GPU CUDA::cudart)
    target_compile_definitions(PIC PRIVATE USE_CUDA)
endif()
endif ()

if(ZLIB_FOUND)
    target_link_libraries(PIC PRIVATE ZLIB::ZLIB)
    target_compile_definitions(PIC PRIVATE HAVE_ZLIB)
endif()

if(PICM_FILESYSTEM_LIB)
    target_link_libraries(PIC PRIVATE ${PICM_FILESYSTEM_LIB})
endif()

if(OpenMP_CXX_FOUND)
    target_link_libraries(PIC PRIVATE OpenMP::OpenMP_CXX)
    target_compile_definitions(PIC PRIVATE USE_OPENMP)
endif()

target_compile_definitions(PIC PRIVATE
    $<IF:$<BOOL:${USE_FLOAT_PRECISION}>,USE_FLOAT,USE_DOUBLE>)

target_compile_definitions(PIC PRIVATE
        $<IF:$<BOOL:${USE_FAST_ALG}>,USE_SPEED,USE_PRECISION>)

target_compile_options(PIC PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:MSVC>:/W4>)

target_compile_options(PIC PRIVATE
    $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:-O3 -march=native>
    $<$<AND:$<CONFIG:Debug>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:-O0 -g>)
