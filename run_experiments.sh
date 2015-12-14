#!/bin/bash

#Let's run the MPI version
for file in x.mm_2D*
do
    if [[ -f $file ]]; then
        PROC_NAME=$file
        cmd="qsub -lnodes=8:ppn=2,mem=9050mb -d . -N $PROC_NAME -e logs -o logs -v EXE=$PROC_NAME,DATA=data/2_repititions run_mm_mpi_pg.sh"
        echo $cmd
        $cmd
    fi
done
