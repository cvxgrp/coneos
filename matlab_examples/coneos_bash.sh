#!/bin/bash
FILES=DIMACS/*
for f in $FILES
do
    foo=${f#*/}
    echo "$foo"
	if [ -f DIMACS_results/${foo}.output ]
	then
		continue
	else
		echo "running test $foo"
        ../coneOSsparse/bin/demo_direct $f > DIMACS_results/$foo.output
	fi
done
