#!/usr/bin/env bash

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
    for i in {1..36}
    do
        (
            cd 6-structure_refinements/output/${datatype}${i}_refined
            phenix.refine --overwrite ../../input/starting_coords.cif ../../input/starting_maps/${datatype}${i}_merged.mtz ../../input/params.eff
        ) &
        pwait $MAX_PROCESSES
    done
done

wait
