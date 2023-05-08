FCOMP = gfortran

FC = -c
FFLAGS = -Wall -Wextra -Wconversion -Wno-maybe-uninitialized -pedantic
FOPT = -O3 -march=native -funroll-loops -flto
FOMP = -fopenmp

EXECUTABLE = run


F90_FILES := $(wildcard src/*.f90)
OBJECTS := $(patsubst 	src/%.f90, src/%.o, $(F90_FILES))

$(EXECUTABLE): $(OBJECTS)
	$(FCOMP) $(FFLAGS) $(FOMP) $(FOPT) $^ -o $@

%.o : %.f90
	$(FCOMP) $(FC) $(FFLAGS) $(FOPT) $(FOMP) $< -o $@

clean:
	rm -rf src/*.o $(EXECUTABLE)