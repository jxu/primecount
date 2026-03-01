CXX 	 := g++
CXXFLAGS := -Wall -Wextra -std=c++17 -fopenmp
SRC_DIR  := src
TEST_DIR := tests
BUILD    := build

# Build specific flags
DEBUG_FLAGS   := -Og -g -D_GLIBCXX_DEBUG -fsanitize=signed-integer-overflow
RELEASE_FLAGS := -O3

# Source files
SRC  := $(wildcard $(SRC_DIR)/*.cpp)
TEST := $(wildcard $(TEST_DIR)/*.cpp)

# Output binaries
BIN_DEBUG   := $(BUILD)/primecount_debug
BIN_RELEASE := $(BUILD)/primecount
BIN_TESTS   := $(BUILD)/tests

# default target
all: release debug

# Release
release: $(BIN_RELEASE)

$(BIN_RELEASE): $(SRC)
	mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) $(RELEASE_FLAGS) $^ -o $@

# Debug
debug: $(BIN_DEBUG)

$(BIN_DEBUG): $(SRC)
	mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) $(DEBUG_FLAGS) $^ -o $@

# Tests (debug mode)
tests: $(BIN_TESTS)

$(BIN_TESTS): $(SRC) $(TEST)
	mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) $(DEBUG_FLAGS) $^ -o $@

clean:
	rm -rf $(BUILD)

.PHONY: all release debug tests clean
