#!/bin/bash
for i in `seq 1 100`;
do
	qsub jobname$i.sh
	sleep 1s
done
