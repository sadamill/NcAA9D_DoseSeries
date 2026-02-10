#!/usr/bin/env bash

# directory setup
mkdir 5-data_processing/input

for i in {1..36}
do
	mkdir -p $(printf "5-data_processing/output/wedge%02d" "$i")
	mkdir -p $(printf "5-data_processing/output/pseudohelix%02d/subwedges" "$i")
done

# copy input files from previous step
cp 4-sampling/output/samples.csv 5-data_processing/input/
