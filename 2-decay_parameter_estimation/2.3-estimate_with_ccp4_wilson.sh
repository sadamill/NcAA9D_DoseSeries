#!/usr/bin/env bash

for start_frame in {1..25}
do
    (
		cd $(printf "2-decay_parameter_estimation/output/dials_output/start_frame_%02d" "$start_frame")
		wilson hklin integrated.mtz <<-EOF > ccp4_wilson.txt
			nresidues 892
			resolution 68.72 1.1
			labin IP=I SIGIP=SIGI
			title wilson.log
		EOF
	)
done

grep "Least squares" start_frame_*/ccp4_wilson.txt | sed -r "s/(^.*B  =  |        SCALE.*$)//g" > 2-decay_parameter_estimation/output/dials_output/ccp4_wilson_b.csv
grep "Least squares" start_frame_*/ccp4_wilson.txt | sed -r "s/(^.*SCALE  =  )//g" > 2-decay_parameter_estimation/output/dials_output/ccp4_wilson_scale.csv

cp 2-decay_parameter_estimation/output/dials_output/*.csv 2-decay_parameter_estimation/input/r_input
