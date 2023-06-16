# Makefile for Project/tests

# Compiler settings
FC = h5fc
FLAGS = -cpp -Wall -Wextra

# Directories
SRC_DIR = ../src
OBJECTS_DIR = objects

# Source files
SRC_FILES = $(wildcard $(SRC_DIR)/*.f90)
TEST_FILES = $(wildcard *.f90)

# Object files
SRC_OBJ_FILES = $(filter-out $(OBJECTS_DIR)/main.o, $(patsubst $(SRC_DIR)/%.f90, $(OBJECTS_DIR)/%.o, $(SRC_FILES)) )
TEST_OBJ_FILES = $(patsubst %.f90,$(OBJECTS_DIR)/%.o,$(TEST_FILES))

# Executable
EXECUTABLE = tests

# Default target
all: $(EXECUTABLE)

# Compile source files
$(OBJECTS_DIR)/%.o: $(SRC_DIR)/%.f90
	@mkdir -p $(OBJECTS_DIR)
	$(FC) $(FLAGS) -c $< -o $@

# Compile test files
$(OBJECTS_DIR)/%.o: %.f90
	@mkdir -p $(OBJECTS_DIR)
	$(FC) $(FLAGS) -c $< -o $@

# Link object files and create executable
$(EXECUTABLE): $(SRC_OBJ_FILES) $(TEST_OBJ_FILES)
	$(FC) $(FLAGS) $^ -o $@

# Clean generated files
clean:
	rm -rf $(OBJECTS_DIR)
	rm -f $(EXECUTABLE)