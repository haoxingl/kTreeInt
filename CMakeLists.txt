cmake_minimum_required(VERSION 3.10)

project(IntegerKTree VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set output directory for executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Add include directories
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR})  # For config.h

# Add source files
set(SOURCES
    src/config.cpp
    src/utilities.cpp
    src/k_sum.cpp
    src/observer.cpp
    src/k_tree.cpp
    src/evaluation.cpp
    src/EvaluationLogger.cpp
    src/ResourceMonitor.cpp
)

# Function to create executable
function(add_project_executable TARGET_NAME)
    add_executable(${TARGET_NAME} ${ARGN} ${SOURCES})
    target_compile_options(${TARGET_NAME} PRIVATE -fdiagnostics-color=always)
    target_compile_options(${TARGET_NAME} PRIVATE $<$<CONFIG:Debug>:-g>)
    target_compile_options(${TARGET_NAME} PRIVATE $<$<CONFIG:Release>:-O3>)
    target_link_libraries(${TARGET_NAME} gmpxx gmp)
endfunction()

# Add main executable
add_project_executable(main ${CMAKE_SOURCE_DIR}/main.cpp)

# Add a custom target for building with O3 optimization
add_custom_target(optimized
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} --build . --target main
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Building main with O3 optimization"
)