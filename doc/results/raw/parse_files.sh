#!/bin/bash

#Let's run the MPI version
for file in x.mm_*
do
    if [[ -f $file ]]; then
        grep "T |" $file > ../$file
    fi
done
