#!/bin/bash
set -ue

echo 'Results prefix: ${1}'
trap "echo Exited!; exit;" SIGINT SIGTERM

for i in $(seq -f "%010g" 1 10000 100000); do 
    jetstream run test_templates/null.jst --int:tasks ${i} --log-level DEBUG --autosave 0 -b > /dev/null &
    JETSTREAM_PID=$!
    psrecord ${JETSTREAM_PID} --interval 1 --log ${1}_${i}.txt || break
done
