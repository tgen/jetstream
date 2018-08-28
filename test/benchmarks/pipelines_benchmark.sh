#!/bin/bash

TASKS=$1
BACKEND=$2
PROJ=$(mktemp -d)

cd ${PROJ}

jetstream project init
jetstream pipelines test/stress_true.yaml --backend ${BACKEND} --no-autosave --int:tasks ${TASKS}

