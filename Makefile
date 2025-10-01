CXX = g++
CXXFLAGS = -Wall -Wextra -fopenmp

# Target-specific variable values
release: CXXFLAGS += -O3
release: primecount

.PHONY: primecount tests

debug: CXXFLAGS += -g -D_GLIBCXX_DEBUG -DDEBUG -fsanitize=signed-integer-overflow 
debug: primecount tests

primecount: main.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

tests: tests.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<
