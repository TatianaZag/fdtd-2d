CC=./dvm c
CFLAGS=-DMEDIUM_DATASET 
all: fdtd-2d-dvmh
	./dvm run 2 2 fdtd-2d-dvmh
fdtd-2d-dvmh: fdtd-2d-dvmh.c
	$(CC) $(CFLAGS) fdtd-2d-dvmh.c
clean:
	rm -f ./fdtd-2d-dvmh ex ey hz
