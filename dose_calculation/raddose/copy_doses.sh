#!/bin/bash

for i in {1..38}
do
    cp ../../../raddose/output/ddwd/wedge${i}/output-DWDs.csv ddwd/wedge_${i}.csv
	cp ../../../raddose/output/fwd/wedge${i}/output-DWDs.csv fwd/wedge_${i}.csv
done