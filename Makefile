CXX = g++
CXXFLAGS = -Wall -Wextra -fopenmp

all: release debug

# Target-specific variable values
release: CXXFLAGS += -O3
release: primecount tests

debug: CXXFLAGS += -Og -g -D_GLIBCXX_DEBUG -fsanitize=signed-integer-overflow
debug: primecount tests

primecount: src/main.cpp src/primecount.hpp src/fenwick_tree.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^

tests: src/tests.cpp src/primecount.hpp src/fenwick_tree.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^

.PHONY: primecount tests
