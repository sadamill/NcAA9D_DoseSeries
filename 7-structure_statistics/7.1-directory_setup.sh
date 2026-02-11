#!/usr/bin/env bash

# copy output files from previous step
for i in {1..36}
do
	cp $(printf "5-data_processing/output/pseudohelix%02d/unmerged.mtz" "$i") 7-structure_statistics/input/unmerged_data/pseudohelix${i}_unmerged.mtz
	cp $(printf "5-data_processing/output/wedge%02d/unmerged.mtz" "$i") 7-structure_statistics/input/unmerged_data/wedge${i}_unmerged.mtz
	cp $(printf "5-data_processing/output/pseudohelix%02d/merged.mtz" "$i") 7-structure_statistics/input/merged_data/pseudohelix${i}_merged.mtz
	cp $(printf "5-data_processing/output/wedge%02d/merged.mtz" "$i") 7-structure_statistics/input/merged_data/wedge${i}_merged.mtz
	cp 6-structure_refinements/output/pseudohelix${i}_refined/Phenix-refine_001.pdb 7-structure_statistics/input/pdbs/pseudohelix${i}.pdb
	cp 6-structure_refinements/output/wedge${i}_refined/Phenix-refine_001.pdb 7-structure_statistics/input/pdbs/wedge${i}.pdb
done