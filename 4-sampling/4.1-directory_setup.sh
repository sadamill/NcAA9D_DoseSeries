#!/usr/bin/env bash

# copy input files from previous step
for i in {1..36}
do
	cp 3-ddwd_calculation/output/r_output/* 4-sampling/input/r_input/
done