#!/usr/bin/env python3
import sys
import jetstream


def main(workflow_path):
    wf = jetstream.workflows.load_workflow(workflow_path)
    tasks = list(wf.tasks(data=True))
    tasks.sort(key=lambda t: t[1]['datetime_start'])

    for task_id, task in tasks:
        method = task.get('methods')
        if method:
            print(method)


if __name__ == "__main__":
    main(sys.argv[1])
