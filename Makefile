CXX = g++
CXXFLAGS = -Wall -Wextra -fopenmp

all: release debug

# Target-specific variable values
release: CXXFLAGS += -O3
release: primecount run_tests

debug: CXXFLAGS += -Og -g -D_GLIBCXX_DEBUG -fsanitize=signed-integer-overflow
debug: primecount run_tests

primecount: src/main.cpp src/primecount.hpp src/fenwick_tree.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^

run_tests: tests/tests.cpp src/primecount.hpp src/fenwick_tree.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^

.PHONY: primecount run_tests
