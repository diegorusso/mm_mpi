#!/bin/bash

# That's the serial version
cmd="./x.mm_serial -f example_data -d"
$cmd

#Let's run the MPI version
for file in x.mm_2D*
do
    if [[ -f $file ]]; then
        cmd="mpirun -np 4 ./$file -f example_data -d"
        $cmd
    fi
done
