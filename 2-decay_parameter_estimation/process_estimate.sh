#!/bin/bash

# directory setup
for i in {1..36}
do
	mkdir -p $(printf "2-decay_parameter_estimation/output/wedge%02d" "$i")
done

for i in {1..25}
do
	mkdir -p $(printf "2-decay_parameter_estimation/output/start_frame_%02d/subwedges" "$i")
done

# import and bin images
for i in {2..37}
do
    (
		new_wedge_num=$((i-1))
	    cd $(printf "2-decay_parameter_estimation/output/wedge%02d" "$new_wedge_num")
	    dials.import template=$(printf "../../../0-diffraction_images/wedge%02d/HR00498_Pn5_wedge_%02d.####.gz" "$i" "$i")
	    dials.find_spots imported.expt
	)
done

for start_frame in {1..25}
do
    (
		cd $(printf "2-decay_parameter_estimation/output/start_frame_%02d" "$start_frame")
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
		dials.export integrated.expt refined.refl
	)
done
