CXX = g++
CXXFLAGS = -std=c++17 -I./include -Wall -Wextra -Wpedantic -Werror -O3 -g

LDFLAGS = -pthread

SRC_DIR = ./src
INCLUDE_DIR = ./include
OBJ_DIR = ./obj

# List of source files
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
# Corresponding object files
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))

# Target executable
TARGET = readon

# Create obj directory if it doesn't exist
$(shell mkdir -p $(OBJ_DIR))

# Build the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)


$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up
clean:
	rm -f $(OBJS) $(TARGET)
