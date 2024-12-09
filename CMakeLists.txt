# Minimum version is CMake 3.26
cmake_minimum_required(VERSION 3.26)

# Export compile commands for the language server
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Project instantiation
project(nbody VERSION 0.2.0.8)

# Raylib config to OpenGL 4.3
set(GRAPHICS GRAPHICS_API_OPENGL_43)

# Add Raylib 5.5
include(FetchContent)
FetchContent_Declare(
    raylib 
    URL https://github.com/raysan5/raylib/archive/refs/tags/5.5.tar.gz)
FetchContent_MakeAvailable(raylib)

# NOTE: This should've been fixed in Raylib 5.5
# Raylib 5.0 is kind of sensitive to this but the upcoming release should
# fix this.
# target_compile_options(raylib PRIVATE -Wno-incompatible-pointer-types)

# Configure the project header
configure_file(include/configuration.txt
    ${PROJECT_SOURCE_DIR}/include/configuration.hpp)

# Build our main executable
add_executable(${PROJECT_NAME}
    src/particle.cpp
    src/main.cpp)

# Use C++20 on target too
set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD_REQUIRED TRUE)
set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 20)

# Include headers here
target_include_directories(${PROJECT_NAME} PRIVATE include)

# Finally link
target_link_libraries(${PROJECT_NAME} raylib)
