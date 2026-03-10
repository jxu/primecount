# Project info
cmake_minimum_required(VERSION 3.16)
project(Primecount LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# CMake has built-in support for OpenMP support in a compiler
# https://cmake.org/cmake/help/latest/module/FindOpenMP.html
find_package(OpenMP REQUIRED)

# Define main functionality as a library (headers optional here)
add_library(primecount_lib
        src/primecount.cpp
        src/phi_block.cpp
)

# Add header files in include path
target_include_directories(primecount_lib PUBLIC src)

# Different builds: Release and Debug
# Debug uses _GLIBCXX_DEBUG which changes all standard class templates
# https://gcc.gnu.org/onlinedocs/libstdc++/manual/debug_mode_using.html
# Debug mode must be applied for all compilation and linking
target_compile_options(primecount_lib PRIVATE
        $<$<CONFIG:Release>:-O3>
        $<$<CONFIG:Debug>:-Og -g -D_GLIBCXX_DEBUG>
)

# Link options: must link in OpenMP
target_link_libraries(primecount_lib PRIVATE OpenMP::OpenMP_CXX)
target_link_options(primecount_lib PRIVATE
        $<$<CONFIG:Debug>:-D_GLIBCXX_DEBUG>
)

# Main executable: links in my library
add_executable(primecount src/main.cpp)
target_link_libraries(primecount PRIVATE primecount_lib)

# Enable CTest (see tests/CMakeLists.txt)
enable_testing()
add_subdirectory(tests)