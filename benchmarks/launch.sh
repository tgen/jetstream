#!/bin/bash
set -eu

for i in $(seq $1 $2 $3); do
  srun -J bm1_${i} -p compute-16core --exclusive bash pipelines_benchmark.sh ${i}
done
	 
