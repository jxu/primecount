CXX = g++
CXXFLAGS = -Wall -Wextra -fopenmp

# Target-specific variable values
release: CXXFLAGS += -O3
release: primecount

debug: CXXFLAGS += -g -D_GLIBCXX_DEBUG -DDEBUG -fsanitize=signed-integer-overflow 
debug: primecount

primecount: primecount.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

