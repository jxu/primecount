CXX = g++
CXXFLAGS = -Wall -Wextra

# Target-specific variable values
release: CXXFLAGS += -O3
release: primecount

debug: CXXFLAGS += -g -D_GLIBCXX_DEBUG -DDEBUG -fsanitize=signed-integer-overflow 
debug: primecount
debug: tests

tests: tests.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

primecount: main.cpp primecount.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<

