#!/bin/bash

for i in {1..38}
do
    cp ../../../../raddose/output/wedges/wedge${i}/output-DWDs.csv wedge_${i}.csv
done