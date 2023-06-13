CCF = -DHDF
FC = -c
FFLAGS = -Wall -Wextra -Wconversion -Wno-maybe-uninitialized -pedantic -ffree-line-length-none
FOPT = -O3 -march=native -funroll-loops -flto
FOMP = -fopenmp

EXECUTABLE = run

F90_FILES := $(wildcard src/*.f90)
SPECIFIC_SOURCE := src/hdf5_utils.f90

ifdef NOHDF5
	FCOMP = gfortran
	F90_FILES := $(filter-out $(SPECIFIC_SOURCE), $(F90_FILES))
	CCF =
else
	FCOMP ?= h5fc
endif

OBJECTS := $(patsubst   src/%.f90, src/%.o, $(F90_FILES))

$(EXECUTABLE): $(OBJECTS)
		$(FCOMP) -cpp $(CCF) $(FFLAGS) $(FOMP) $(FOPT) $^ -o $@

%.o : %.f90
		$(FCOMP) -cpp $(CCF) $(FC) $(FFLAGS) $(FOPT) $(FOMP) $< -o $@

clean:
		rm -rf src/*.o $(EXECUTABLE)