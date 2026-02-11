#!/usr/bin/env bash

# make destination folders for input and output files
mkdir 7-structure_analysis/input/unmerged_data

# copy output files from previous step
for i in {1..36}
do
	cp $(printf "5-data_processing/output/pseudohelix${i}/merged.mtz") 6-structure_refinements/input/starting_maps/pseudohelix${i}_merged.mtz
	cp $(printf "5-data_processing/output/wedge${i}/merged.mtz") 6-structure_refinements/input/starting_maps/wedge${i}_merged.mtz
done