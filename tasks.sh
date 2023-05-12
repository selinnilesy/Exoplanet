#!/bin/bash
files=(mockcurve*.fits)
for i in  0 1 2 3 4 5
do
	for f in "${files[@]}";
	do
	   export OMP_NUM_THREADS=$((2 **i)) &&
	   echo "- Processing file: $f" &&
		./ompBLS $f > "./logs/wo-mem/log-$f-$((2 **i)).txt"
	done
done