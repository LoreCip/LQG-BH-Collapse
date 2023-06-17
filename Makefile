CCF = -DHDF
FC = -c
FFLAGS = -Wall -Wextra -Wconversion -Wno-maybe-uninitialized -pedantic -ffree-line-length-none
FOPT = -O3 -march=native -funroll-loops -flto
FOMP = -fopenmp

EXECUTABLE = run

OBJECTS_DIR = src/objects

F90_FILES := $(wildcard src/*.f90)
SPECIFIC_SOURCE := src/hdf5_utils.f90

ifdef NOHDF5
	FCOMP = gfortran
	F90_FILES := $(filter-out $(SPECIFIC_SOURCE), $(F90_FILES))
	CCF =
else
	FCOMP ?= h5fc
endif

OBJECTS := $(patsubst   src/%.f90, $(OBJECTS_DIR)/%.o, $(F90_FILES))

$(EXECUTABLE): $(OBJECTS)
		$(FCOMP) -cpp $(CCF) $(FFLAGS) $(FOMP) $(FOPT) $^ -o $@

$(OBJECTS_DIR)/%.o : src/%.f90
		@mkdir -p $(OBJECTS_DIR)
		$(FCOMP) -cpp $(CCF) $(FC) $(FFLAGS) $(FOPT) $(FOMP) $< -o $@

clean:
		rm -rf $(OBJECTS_DIR) 
		rm -rf $(EXECUTABLE)