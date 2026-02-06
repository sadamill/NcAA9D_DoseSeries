#!/bin/bash

pwait() {
    while [[ $(jobs | wc -l) -ge $1 ]]
    do
        wait -n
    done
}

MAX_PROCESSES=12
datatypes=("pseudohelix" "wedge")

for datatype in ${datatypes[@]}
do
    for i in $(seq 1 36)
    do
        (
            cd ./${datatype}${i}/refine_1/
            phenix.refine --overwrite ../../starting_coords.pdb ../starting_map/${datatype}${i}_merged.mtz ../../params.eff
        ) &
        pwait $MAX_PROCESSES
    done
done

wait
