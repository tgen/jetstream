#!/bin/bash
set -ue

# Two required args: a results file prefix name, and path to a template file
# the template should handle one variable "tasks" an integer which defines
# the number of tasks to add. Tests will run for several values.

echo "Results prefix: ${1}"
echo "Template path: ${2}"
trap "echo Exited!; exit;" SIGINT SIGTERM

for i in $(seq -f "%010g" 100001 10000 200000); do 
    jetstream run --build-only ${2} -- --int:tasks ${i} > /dev/null &
    JETSTREAM_PID=$!
    psrecord ${JETSTREAM_PID} --interval 1 --log ${1}_${i}.txt || break
done
