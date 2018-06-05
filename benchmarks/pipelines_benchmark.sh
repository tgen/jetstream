#!/bin/bash

TASKS=$1
PROJ=$(mktemp -d)

cd ${PROJ}

jetstream project init
jetstream pipelines test/stress_true.yaml --int:tasks ${TASKS}

