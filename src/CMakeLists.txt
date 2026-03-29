file(GLOB SOURCES "main.cpp" "core/*.cpp" "solvers/PIC/*.cpp" "solvers/SemiLagrangian/*.cpp")
set(CMAKE_NINJA_FORCE_RESPONSE_FILE "ON" CACHE BOOL "Force Ninja to use response files.")

add_executable(PIC ${SOURCES})
set_target_properties(PIC PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                     "${CMAKE_BINARY_DIR}/bin")

# mandatory library ( downloaded if not available )
target_link_libraries(PIC PRIVATE nlohmann_json::nlohmann_json)

# for compression of VTK files
if(ZLIB_FOUND)
  target_link_libraries(PIC PRIVATE ZLIB::ZLIB)
  target_compile_definitions(PIC PRIVATE HAVE_ZLIB)
endif()
# same for openMP
if(OpenMP_CXX_FOUND)
  target_link_libraries(PIC PRIVATE OpenMP::OpenMP_CXX)
  target_compile_definitions(PIC PRIVATE USE_OPENMP)
endif()

if(USE_FP16)
    if(USE_FLOAT_PRECISION)
        message(FATAL_ERROR "Cannot use USE_FP16 and USE_FLOAT_PRECISION together")
    endif()
    target_compile_definitions(PIC PRIVATE USE_HALF)
elseif(USE_FLOAT_PRECISION)
    target_compile_definitions(PIC PRIVATE USE_FLOAT)
else()
    target_compile_definitions(PIC PRIVATE USE_DOUBLE)
endif()
# vebose build
#target_compile_options(
  #PIC PRIVATE $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-O3 -march=native -Wall -Wextra -Wpedantic>
              #$<$<CXX_COMPILER_ID:MSVC>:/W4>)
# Warnings
target_compile_options(PIC PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>
        $<$<CXX_COMPILER_ID:MSVC>:/W4>
)

# Release optimizations
target_compile_options(PIC PRIVATE
        $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:-O3 -march=native>
)

# Debug sanitizers
target_compile_options(PIC PRIVATE
        $<$<AND:$<CONFIG:Debug>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:
        -O0 -g -fsanitize=address,undefined,leak -fno-omit-frame-pointer
        >
)
target_link_options(PIC PRIVATE
        $<$<AND:$<CONFIG:Debug>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:
        -fsanitize=address,undefined,leak
        >
)