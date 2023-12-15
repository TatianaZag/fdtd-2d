CC=pgcc
CFLAGS=-acc -ta=nvidia -Wall -Wextra -DLARGE_DATASET
all: fdtd-2d-oacc
	./fdtd-2d-oacc
fdtd-2d-oacc: fdtd-2d-oacc.c
	$(CC) $(CFLAGS) fdtd-2d-oacc.c -o fdtd-2d-oacc
clean:
	rm -f ./fdtd-2d-oacc ex ey hz
