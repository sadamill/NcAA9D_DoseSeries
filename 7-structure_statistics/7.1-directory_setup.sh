#!/usr/bin/env bash

# copy output files from previous step
for i in {1..36}
do
	cp $(printf "5-data_processing/output/pseudohelix${i}/unmerged.mtz") 7-structure_statistics/input/unmerged_data/pseudohelix${i}_unmerged.mtz
	cp $(printf "5-data_processing/output/wedge${i}/unmerged.mtz") 7-structure_statistics/input/unmerged_data/wedge${i}_unmerged.mtz
	cp $(printf "5-data_processing/output/pseudohelix${i}/merged.mtz") 7-structure_statistics/input/merged_data/pseudohelix${i}_merged.mtz
	cp $(printf "5-data_processing/output/wedge${i}/merged.mtz") 7-structure_statistics/input/merged_data/wedge${i}_merged.mtz
	cp $(printf "6-structure_refinements/output/pseudohelix${i}_refined/Phenix-refine_001.pdb") 7-structure_statistics/input/pdbs/pseudohelix${i}.pdb
	cp $(printf "6-structure_refinements/output/wedge${i}_refined/Phenix-refine_001.pdb") 7-structure_statistics/input/pdbs/wedge${i}.pdb
done