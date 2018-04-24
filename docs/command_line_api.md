
## Commands

The installation process will make two commands available:

  - `jetstream` Provides access to utility/core functions for building projects and workflows.

  - `jetstream_pipelines` Runs the built-in analysis pipelines.

View the help with `-h/--help` to get started. If you receive an error that acommand was not found, the Python packages bin location is probably not added to your `$PATH`. [Refer to this post][path_help] for more help.

## Scripts

In addition to the two commands mentioned above, individual analysis modules are available to run as scripts and will be added to your `$PATH`. For example,`make_exome_refkit.py -h` can be launched anywhere with no additional parameters needed. 


## Python API

`>>> import jetstream`

TODO: Build and link module docs, maybe readthedocs and github.io?