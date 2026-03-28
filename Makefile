# Compiler
CXX := g++

# Directories
SRC_DIR := src
BUILD_DIR := build
EIGEN_DIR := eigen-5.0.1

# Executable name
TARGET := kpm.out

# Compiler flags
CXXFLAGS := -std=c++23 -O3 -mavx -Wall -Wextra -Wpedantic -DNDEBUG -I$(SRC_DIR) -I$(EIGEN_DIR)

# Linker flags
LDFLAGS :=

# Source and object files
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRCS))

# Default target
all: $(TARGET)

# Link executable
$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

# Compile source files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Create build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Clean build files
clean:
	rm -r $(BUILD_DIR) $(TARGET)

# Rebuild from scratch
rebuild: clean all

# Run the program
run: $(TARGET)
	./$(TARGET)

.PHONY: all clean rebuild run
