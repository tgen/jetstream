#!/bin/bash

module load python

for f in $(seq 0 50000 1000000); do 
  srun -c1 -n1 -J wf_${f} python3 workflows_benchmark2.py ${f} > test_wf_${f}.log 2>&1 &
  srun -c1 -n1 -J pl_${f} bash pipelines_benchmark.sh ${f} local > test_pl_${f}.log 2>&1 &
done

wait

