#!/bin/bash

for id in 1 8 16 32 256
do
make clean;
CFLAGS=-DNUM_TASKS=$id make
done