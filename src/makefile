CC = mpicc
CFLAGS = -O2
DEPS = global.h
OBJ = poisson.o readGrid.o initializationSteady.o cellNeighbours.o faceDistance.o solver.o nodes2Cells.o writingVTK.o

%.o: %.C $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

poissonSolver: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f *.o poissonSolver
