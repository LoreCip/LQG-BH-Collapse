FCOMP = gfortran

FC = -c
FFLAGS = -Wall -Wextra -Wconversion -pedantic
FOPT = -O3
FOMP = -fopenmp

EXECUTABLE = run

FDEBUG = #-g -fcheck=bounds

F90_FILES := $(wildcard src/*.f90)
OBJECTS := $(patsubst 	src/%.f90, src/%.o, $(F90_FILES))

$(EXECUTABLE): $(OBJECTS) $(MODFILES)
	$(FCOMP) $(FDEBUG) $(FFLAGS) $(FOMP) $(FOPT) $^ -o $@

%.o %.mod: %.f90
	$(FCOMP) $(FC) $(FDEBUG) $(FFLAGS) $(FOPT) $(FOMP) $< -o $@

clean:
	rm -rf src/*.o src/*.mod $(EXECUTABLE)