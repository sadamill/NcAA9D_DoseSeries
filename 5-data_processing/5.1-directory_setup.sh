#!/usr/bin/env bash

# make empty input directory
mkdir 5-data_processing/input


# copy input files from previous step
cp 4-sampling/output/r_output/samples.csv 5-data_processing/input/
