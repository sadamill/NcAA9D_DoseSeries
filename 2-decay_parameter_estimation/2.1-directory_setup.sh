#!/usr/bin/env bash

# directory setup
mkdir 2-decay_parameter_estimation/output/r_output

for i in {1..36}
do
	mkdir -p $(printf "2-decay_parameter_estimation/output/dials_output/wedge%02d" "$i")
done

for i in {1..25}
do
	mkdir -p $(printf "2-decay_parameter_estimation/output/dials_output/start_frame_%02d/subwedges" "$i")
done

# copy input files from previous step
cp 1-fwd_calculation/output/r_output/fwds.csv 2-decay_parameter_estimation/input/r_input/
