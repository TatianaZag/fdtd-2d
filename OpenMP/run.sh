#!/bin/bash

#`
for dat in -DSMALL_DATASET -DMEDIUM_DATASET -DLARGE_DATASET
do
    for num in 1 2 4 8 12 20
    do
    echo
    echo "-------< Dataset: $dat count tasks: $num >-------"
    echo
    make clean;
        for ((i = 0; i < 10; i++))
        do #for i
        CFLAGS="$dat -DNUM_TASKS=$num" make
        done #for i
    done #for num
done #for dat
