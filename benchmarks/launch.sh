#!/bin/bash
set -eu

for i in $(seq $1 $2 $3); do
  sbatch -J bm1_${i} -p compute-16core --exclusive benchmark.sh ${i}
done
	 
