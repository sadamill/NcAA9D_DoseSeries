#!/usr/bin/env bash

# directory setup
for i in {1..38}
do
	mkdir -p 3-ddwd_calculation/output/raddose_output/wedge$i
done

mkdir  3-ddwd_calculation/output/r_output

# copy input files from previous step
cp 2-decay_parameter_estimation/output/r_output/decay_params.csv 3-ddwd_calculation/input/raddose_input/
cp 1-fwd_calculation/output/r_output/fwds.csv 3-ddwd_calculation/input/r_input