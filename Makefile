CFLAGS= -Wall -Wextra -DMEDIUM_DATASET -fopenmp
all: fdtd-2d-omp fdtd-2d-tasks
	./fdtd-2d-omp
	./fdtd-2d-tasks
fdtd-2d-omp: fdtd-2d-omp-1.c
fdtd-2d-tasks: fdtd-2d-tasks.c
clean:
	rm -f ./fdtd-2d-omp ./fdtd-2d-tasks