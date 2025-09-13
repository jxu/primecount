CXX = g++
CXXFLAGS = -Wall -Wextra

testsrc = test_fenwick_tree test_phi_block

# Target-specific variable values
release: CXXFLAGS += -O3
release: primecount

debug: CXXFLAGS += -g -D_GLIBCXX_DEBUG -DDEBUG -fsanitize=signed-integer-overflow 
debug: primecount
debug: $(testsrc)

$(testsrc): %: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

primecount: main.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^

