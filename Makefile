CXX 	 := g++
CXXFLAGS := -Wall -Wextra -std=c++17 -fopenmp
SRC_DIR		:= src
TEST_DIR	:= tests
BUILD_DIR	:= build

# Build specific flags
DEBUG_FLAGS   := -Og -g -D_GLIBCXX_DEBUG -fsanitize=signed-integer-overflow
RELEASE_FLAGS := -O3

# Source files
SRC  := src/phi_block.cpp src/primecount.cpp
TEST := $(wildcard $(TEST_DIR)/*.cpp)

# Output binaries
BIN_DEBUG   := $(BUILD_DIR)/primecount_debug
BIN_RELEASE := $(BUILD_DIR)/primecount
BIN_TESTS   := $(BUILD_DIR)/tests

# default target
all: release debug tests

# Release
release: $(BIN_RELEASE)

$(BIN_RELEASE): $(SRC) $(SRC_DIR)/main.cpp
	mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(RELEASE_FLAGS) $^ -o $@

# Debug
debug: $(BIN_DEBUG)

$(BIN_DEBUG): $(SRC) $(SRC_DIR)/main.cpp
	mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(DEBUG_FLAGS) $^ -o $@

# Tests (debug mode)
tests: $(BIN_TESTS)

$(BIN_TESTS): $(SRC) $(TEST)
	mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(DEBUG_FLAGS) $^ -o $@

clean:
	rm -rf $(BUILD_DIR)

.PHONY: all release debug tests clean
