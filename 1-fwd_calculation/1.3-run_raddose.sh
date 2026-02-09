#!/usr/bin/env bash

# define function for parallel processing
pwait() {
    while [[ $(jobs | wc -l) -ge $1 ]]
    do
        wait -n
    done
}

MAX_PROCESSES=8

# calculate FWDs in parallel using pwait
for i in {1..38}
do
    (
    	cd 1-fwd_calculation/output/raddose_output/wedge$i
    	java -jar /Applications/RADDOSE-3D-master/raddose3d.jar -i ../../../input/raddose_input/wedge${i}.txt
    ) &
    pwait $MAX_PROCESSES
done

wait

# move raddose results to folder for R input
for i in {1..38}
do
	cp 1-fwd_calculation/output/raddose_output/wedge${i}/output-DWDs.csv 1-fwd_calculation/input/r_input/wedge${i}.csv
done
