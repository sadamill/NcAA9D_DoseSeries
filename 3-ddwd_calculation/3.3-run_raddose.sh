#!/usr/bin/env bash

# define function for parallel processing
pwait() {
    while [[ $(jobs | wc -l) -ge $1 ]]
    do
        wait -n
    done
}

MAX_PROCESSES=4

# calculate fwds in parallel using pwait
for i in {1..38}
do
    (
    	cd 3-ddwd_calculation/output/raddose_output/wedge$i
    	java -jar /Applications/RADDOSE-3D-master/raddose3d.jar -i ../../../input/raddose_input/wedge${i}.txt
    ) &
    pwait $MAX_PROCESSES
done

wait

# copy output files to the r input
for i in {1..38};
do
	cp 3-ddwd_calculation/output/raddose_output/wedge$i/output-DWDs.csv 3-ddwd_calculation/input/r_input/wedge${i}.csv
done
