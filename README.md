
# Installation

## Make sure you have Python3 installed

If you already have Python3 and Pip installed, skip to the next step.

This requires Python3. For a lot of reasons, [Homebrew](https://brew.sh) is 
the best way to install Python3. Follow this guide:

http://docs.python-guide.org/en/latest/starting/install3/osx/

## Single Command Install with Pip

# Development notes

_Project_

A directory that contains a `.jetstream` index. Projects also contain a
collection of data that can be analyzed with jetstream. Projects can be
initialized by jetstream, or existing projects can be updated by jetstream.

_Run_

A singular instance of the jetstream application operating on a project. Each
run creates a new record in the project index. Each run is uniquely identified
by its ID, a 26 character string. `jetstream resume` is able to restart a run
that is not complete.

_Index_

The `.jetstream` directory, and its contents, comprise a jetstream index.

_Workflow_

The set of plugins that will be executed during a run. This is organized as a
directed acyclic graph to respect plugin dependencies.

_Config_

Config file, run config, etc.. This is a document describing workflow and
settings to use for a run.


## Slurm remote execution

Should slurm strategy allow for the head node to be a remote machine? This
would enable the Jetstream controller to run anywhere and issue jobs to a
slurm scheduler.

To do this, a local data store would need to be created on a shared filesystem
or a remote datastore would need to be established.

Heres an example of launching a slurm job over ssh:

`cat test_slurm.sh | ssh dback-login1.tgen.org 'bash -l -c "srun bash"'`

This would require the application to be configured with the address and
potentially login credentials for the slurm head node.

# Extra features

There are several general purpose genomics modules included in this package.
Here is an example of processing intervals files with functions included in
`jetstream.formats.intervals`: 

```python
from jetstream.formats.intervals import read_gffv2, to_bed

ints = intervals.read_gffv2('targets.gtf')
ints = ints.filter(lambda i: i['feature'] == 'exon')  # Select only exons

with open('filtered_targets.bed', 'w') as fp:
    fp.write(intervals.to_bed(ints))
```
