add_executable(PIC main.cpp)
# Optional: Set output directory
set_target_properties(PIC PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
