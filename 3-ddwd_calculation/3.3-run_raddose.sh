#!/usr/bin/env bash

# directory setup
for i in {1..38}
do
	mkdir -p 3-ddwd_calculation/output/raddose_output/wedge_$i
done

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
    	cd 3-ddwd_calculation/output/raddose_output/wedge_$i
    	java -jar /Applications/RADDOSE-3D-master/raddose3d.jar -i ../../../input/raddose_input/wedge${i}.txt
    ) &
    pwait $MAX_PROCESSES
done

wait
