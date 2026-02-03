# Optional: Set output directory
set(SOURCES main.cpp core/Grid2D.cpp core/Parameters.cpp core/OutputWriter.cpp)

add_executable(PIC ${SOURCES})
set_target_properties(PIC PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                     "${CMAKE_BINARY_DIR}/bin")

target_link_libraries(PIC PRIVATE Eigen3::Eigen)
# JSON STYLE
target_link_libraries(PIC PRIVATE nlohmann_json::nlohmann_json)
# kokkos ? target_link_libraries(MyTarget PRIVATE Kokkos::kokkos) Compiler
# warnings
target_compile_options(PIC PRIVATE $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra
                                   -Wpedantic> $<$<CXX_COMPILER_ID:MSVC>:/W4>)
