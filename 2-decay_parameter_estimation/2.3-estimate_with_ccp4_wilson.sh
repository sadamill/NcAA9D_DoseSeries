#!/usr/bin/env bash

for start_frame in {1..25}
do
    cd $(printf "start_frame_%02d" "$start_frame")
	wilson hklin integrated.mtz <<-EOF > ccp4_wilson.txt
	nresidues 892
	resolution 68.72 1.1
	labin IP=I SIGIP=SIGI
	title wilson.log
	EOF
    cd ..
done

grep "Least squares" start_frame_*/ccp4_wilson.txt | sed -r "s/(^.*B  =  |        SCALE.*$)//g" > ccp4_wilson_b.csv
grep "Least squares" start_frame_*/ccp4_wilson.txt | sed -r "s/(^.*SCALE  =  )//g" > ccp4_wilson_scale.csv
