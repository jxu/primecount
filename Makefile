CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++20 -fopenmp

# Target-specific variable values
release: CXXFLAGS += -O3
release: primecount

.PHONY: primecount tests

debug: CXXFLAGS += -Og -g -D_GLIBCXX_DEBUG -DDEBUG -fsanitize=signed-integer-overflow
debug: primecount tests

primecount: main.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

tests: tests.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<
