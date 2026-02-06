#!/usr/bin/env bash

# copy input files from previous step
for i in {1..36}
do
	cp 2-decay_parameter_estimation/output/r_output/*___________.csv 3-decay_parameter_estimation/input/raddose_input/
done