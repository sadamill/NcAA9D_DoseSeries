#!/usr/bin/env bash

# make destination folders for input and output files
mkdir 6-structure_refinements/input/starting_maps

for i in {1..36}
do
	mkdir -p 6-structure_refinements/output/wedge${i}_refined 
	mkdir -p 6-structure_refinements/output/pseudohelix${i}_refined
done

# copy output files from previous step
for i in {1..36}
do
	cp $(printf "5-data_processing/output/pseudohelix%02d/merged.mtz" "$i") 6-structure_refinements/input/starting_maps/pseudohelix${i}_merged.mtz
	cp $(printf "5-data_processing/output/wedge%02d/merged.mtz" "$i") 6-structure_refinements/input/starting_maps/wedge${i}_merged.mtz
done