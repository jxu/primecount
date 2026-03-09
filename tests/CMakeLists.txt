# Collect all test source files in this directory
file(GLOB TEST_SOURCES *.cpp)

# Create the test executable
add_executable(primecount_tests ${TEST_SOURCES})

# Link against your library and OpenMP
target_link_libraries(primecount_tests PRIVATE primecount_lib OpenMP::OpenMP_CXX)

# Register with CTest
add_test(NAME primecount_unit_tests COMMAND primecount_tests)