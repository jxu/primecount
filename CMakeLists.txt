# Project info
cmake_minimum_required(VERSION 3.16)
project(Primecount LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# CMake has built-in support for OpenMP support in a compiler
find_package(OpenMP REQUIRED)

# Define main functionality as a library
add_library(primecount_lib
        src/primecount.cpp
        src/phi_block.cpp
)

target_include_directories(primecount_lib PUBLIC src)

# Different builds: Release and Debug
target_compile_options(primecount_lib PRIVATE
        $<$<CONFIG:Release>:-O3>
        $<$<CONFIG:Debug>:-Og -g -D_GLIBCXX_DEBUG -fsanitize=signed-integer-overflow>
)

target_link_options(primecount_lib PRIVATE OpenMP::OpenMP_CXX)

# Main executable
add_executable(primecount src/main.cpp)
target_link_libraries(primecount PRIVATE primecount_lib OpenMP::OpenMP_CXX)

# Tests directory (see tests/CMakeLists.txt)
enable_testing()
add_subdirectory(tests)