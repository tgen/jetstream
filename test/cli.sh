#!/usr/bin/env bash

set -ue

run_cmd(){
    CMD="$@"
    printf "TESTING: ${CMD}\n"
    $@
}

DIR=$(mktemp -d)
cd ${DIR}


run_cmd jetstream -v
run_cmd jetstream -h
run_cmd jetstream init test_project

cd test_project

run_cmd jetstream run <(echo '- cmd: hostname')
run_cmd jetstream pipelines <(echo '- cmd: hostname')
run_cmd jetstream project config
run_cmd jetstream project tasks

rm -r ${DIR}