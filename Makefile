CXX = g++
CXXFLAGS = -Wall -Wextra -fopenmp

# Target-specific variable values
release: CXXFLAGS += -O3
release: primecount tests

debug: CXXFLAGS += -Og -g -D_GLIBCXX_DEBUG -fsanitize=signed-integer-overflow
debug: primecount tests

primecount: main.cpp primecount.hpp fenwick_tree.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^

tests: tests.cpp primecount.hpp fenwick_tree.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^

.PHONY: primecount tests
