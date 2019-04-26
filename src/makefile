CC = gcc
CFLAGS = -lm -fopenmp -O2
CFLAGS2 = -lm -fopenmp -lpthread -O2
DEPS = global.h
OBJ = poisson.o readGrid.o initializationSteady.o cellNeighbours.o faceDistance.o solver.o nodes2Cells.o writingVTK.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

poissonSolver: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS2)

.PHONY: clean

clean:
	rm -f *.o poissonSolver
