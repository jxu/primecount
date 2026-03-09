CXX 	    := g++
CXXFLAGS    := -Wall -Wextra -std=c++17 -fopenmp
LDFLAGS     := -fopenmp
SRC_DIR		:= src
TEST_DIR	:= tests
BUILD_DIR	:= build

# Build mode: Debug or Release specific flags
MODE ?= Release

ifeq ($(MODE),Release)
	CXXFLAGS += -O3
else
	CXXFLAGS += -Og -g -D_GLIBCXX_DEBUG -fsanitize=signed-integer-overflow
	LDFLAGS += -fsanitize=signed-integer-overflow
endif

# Source files
#SRC  := src/phi_block.cpp src/primecount.cpp
SRC := $(wildcard $(SRC_DIR)/*.cpp)

# Output binaries
BIN_TARGET  := $(BUILD_DIR)/primecount
BIN_TESTS   := $(BUILD_DIR)/tests

# Object files go into build/
OBJ := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRC))

# default target
all: $(BIN_TARGET)

$(BIN_TARGET): $(OBJ) $(SRC_DIR)/main.cpp
	$(CXX) $(LDFLAGS) $(OBJ) -o $@

# Compile rule for object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@


# Tests (debug mode)
TEST_SRC := $(wildcard tests/*.cpp)
TEST_OBJ := $(patsubst tests/%.cpp,$(BUILD_DIR)/tests_%.o,$(TEST_SRC))

tests: $(BIN_TESTS)

$(BIN_TESTS): $(TEST_OBJ) $(OBJ)
	$(CXX) $(TEST_OBJ) $(OBJ) -o tests

clean:
	rm -rf $(BUILD_DIR)

.PHONY: all tests clean
