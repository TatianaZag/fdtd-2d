#!/bin/bash

N=10
for dat in -DSMALL_DATASET -DMEDIUM_DATASET -DLARGE_DATASET
do
    for num in 1 2 4 8 12 20
    do
    echo
    echo "-------< Dataset: $dat count tasks: $num >-------"
    echo
    make clean;
        sum=0
        for ((i = 0; i < $N; i++))
        do #for i
        res=$(CFLAGS="$dat -DNUM_TASKS=$num" make | grep -Eo "[0-9]+\.[0-9]+")
        sum=$(echo $res + $sum | bc)
        done #for i
        echo "scale=6; $sum/$N" | bc
    done #for num
done #for dat
