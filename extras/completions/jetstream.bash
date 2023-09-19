#!/bin/bash

_jetstream_compgen_filenames() {
    local cur="$1"

    # Files, excluding directories:
    grep -v -F -f <(compgen -d -P ^ -S '$' -- "$cur") \
        <(compgen -f -P ^ -S '$' -- "$cur") |
        sed -e 's/^\^//' -e 's/\$$/ /'

    # Directories:
    compgen -d -S / -- "$cur"
}

_jetstream() {
    local cur prev opts

    cur=${COMP_WORDS[COMP_CWORD]}
    prev=${COMP_WORDS[COMP_CWORD-1]}
    opts="--mash -m --backend --build -b --render -r --status -s --output -o --verbose -v --logging -l --no-logs"

    if [[ ${prev} == --backend ]]; then
        COMPREPLY=( $(compgen -W "slurm slurm_singularity" -- ${cur}) )
        return 0
    fi

    if [[ ${cur} == -* ]]; then
        COMPREPLY=( $(compgen -W "${opts}" -- ${cur}) )
        return 0
    fi

    case ${COMP_CWORD} in
        1)
            COMPREPLY=($(compgen -W "init run mash build render tasks project pipelines settings" -- ${cur}))
            ;;
        2)
            case ${prev} in
                pipelines)
                    COMPREPLY=($(compgen -W "$(head ${HOME}/jetstream_pipelines/*/pipeline.yaml | grep '^  name:\|^  version:' | awk '{ print $2 }' | xargs -n2 | sed 's/ /@/g')" -- ${cur}))
                    ;;
            esac
            ;;
        *)
            COMPREPLY=($(_jetstream_compgen_filenames "$cur"))
            ;;
    esac
}

complete -o nospace -F _jetstream jetstream
