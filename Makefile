FC = gfortran
FFLAGS = -g -c -O3 -Wall -Wextra -Wconversion -pedantic -fcheck=bounds -fmax-errors=5
EXECUTABLE = run

F90_FILES := $(wildcard src/*.f90)
OBJECTS := $(patsubst 	src/%.f90, src/%.o, $(F90_FILES))

$(EXECUTABLE): $(OBJECTS) $(MODFILES)
	$(FC) $^ -o $@

%.o %.mod: %.f90
	$(FC) $(FFLAGS) $< -o $@

clean:
	rm -rf *.o *.mod $(EXECUTABLE)