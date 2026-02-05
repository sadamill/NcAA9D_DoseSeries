#!/bin/bash

for i in {1..36}
do
    mkdir -p pseudohelix$i/starting_map pseudohelix$i/refine_1
    cp ../data_processing/output/pseudohelix${i}_merged.mtz pseudohelix$i/starting_map
    mkdir -p wedge$i/starting_map wedge$i/refine_1
    cp ../data_processing/output/wedge${i}_merged.mtz wedge$i/starting_map/wedge${i}_merged.mtz
done