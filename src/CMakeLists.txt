file(GLOB SOURCES "main.cpp" "core/*.cpp" "solvers/SemiLagrangian/*.cpp")
set(CMAKE_NINJA_FORCE_RESPONSE_FILE "ON" CACHE BOOL "Force Ninja to use response files.")

add_executable(PIC ${SOURCES})
set_target_properties(PIC PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# not necesarry anymore i suppose
#target_include_directories(PIC PRIVATE /usr/include/paraview)

target_link_libraries(PIC PRIVATE 
    nlohmann_json::nlohmann_json
    OpenMP::OpenMP_CXX
)
# include lib thus
if(ZLIB_FOUND)
  target_link_libraries(PIC PRIVATE ZLIB::ZLIB)
  target_compile_definitions(PIC PRIVATE HAVE_ZLIB)
endif()

if(USE_FLOAT_PRECISION)
    target_compile_definitions(PIC PRIVATE USE_FLOAT)
else()
    target_compile_definitions(PIC PRIVATE USE_DOUBLE)
endif()
target_compile_options(PIC PRIVATE 
    $<$<CXX_COMPILER_ID:GNU,Clang>:-O3 -march=native -Wall -Wextra -Wpedantic> 
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
)
