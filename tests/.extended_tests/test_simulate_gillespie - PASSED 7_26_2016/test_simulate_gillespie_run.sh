#!/bin/sh

for sim_num in 1 2 3; do
	#
	echo "${sim_num}"
	export sim_num
	#
	sbatch --job-name=test_gillespie_${sim_num} \
	test_gillespie.sbatch
	#
	sleep 1 #
done
