#!/usr/bin/env bash

# import and process wedges
for i in {2..37}
do
    (
		new_wedge_num=$((i-1))
	    cd $(printf "5-data_processing/output/wedge%02d" "$new_wedge_num")
	    dials.import template=$(printf "../../../0-diffraction_images/wedge%02d/HR00498_Pn5_wedge_%02d.####.gz" "$i" "$i")
	    dials.find_spots imported.expt
	    dials.index imported.expt strong.refl space_group=P1211
	    dials.refine indexed.expt indexed.refl
	    dials.integrate refined.expt refined.refl
	    dials.scale integrated.expt integrated.refl unmerged_mtz=unmerged.mtz merged_mtz=merged.mtz d_min=1.1
	)
done

# start frames obtained from weighted random sampling
start_frames=(1 2 3 4 5 6 8 9 10 11 13 15 20 23 26 27 32 37 51 76 77 92 101 118 126 149 151 166 167 168 169 171 173 174 175 176)

# iterate through indices of array to process pseudohelices
for i in ${!start_frames[@]}
do
	(
		pseudohelix_number=$((i+1))
	    cd $(printf "5-data_processing/output/pseudohelix%02d" "$pseudohelix_number")
	    start_frame=$((${start_frames[$i]}))
	    end_frame=$((${start_frame}+4))
	    for wedge in {1..36}
	    do
	        dials.slice_sequence $(printf "../wedge%02d/imported.expt" "$wedge") $(printf "../wedge%02d/strong.refl" "$wedge") \
	        					 "image_range=$start_frame $end_frame" \
	        					 reflections_filename=$(printf "subwedges/wedge%02d_slice_${start_frame}-${end_frame}.refl"  "$wedge") \
	        					 experiments_filename=$(printf "subwedges/wedge%02d_slice_${start_frame}-${end_frame}.expt" "$wedge")
	    done
	    dials.combine_experiments subwedges/wedge*_slice_${start_frame}-${end_frame}.expt subwedges/wedge*_slice_${start_frame}-${end_frame}.refl
	    dials.index combined.expt combined.refl space_group=P1211
	    dials.refine indexed.expt indexed.refl
	    dials.integrate refined.expt refined.refl
	    dials.scale integrated.expt integrated.refl unmerged_mtz=unmerged.mtz merged_mtz=merged.mtz d_min=1.1
	)
done
