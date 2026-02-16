#!/usr/bin/env bash

# copy input files over
for i in {1..36}
do
	cp 6-structure_refinements/output/pseudohelix${i}_refined/Phenix_refine_001.pdb 8-structure_analysis/input/structures/pseudohelix${i}.pdb
	cp 6-structure_refinements/output/wedge${i}_refined/Phenix_refine_001.pdb 8-structure_analysis/input/structures/wedge${i}.pdb
done

cp 7-structure_statistics/output/* 8-structure_analysis/input/crystal_statistics/
cp 4-sampling/output/samples.csv 8-structure_analysis/input/
cp 3-ddwd_calculation/output/r_output/ddwds.csv 8-structure_analysis/input