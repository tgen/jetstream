Command-Line
============


    ``jetstream``

    All command-line features can be accessed with ``jetstream <command>``
    Arguments described in this section can be used with any command.

    .. argparse::
        :ref: jetstream.cli.jetstream.arg_parser
        :prog: jetstream


pipelines
----------

    ``jetstream pipelines``


    .. argparse::
        :ref: jetstream.cli.subcommands.pipelines.arg_parser
        :prog: jetstream pipelines
        :nodefault:

--------

project
--------

    ``jetstream project``


    .. argparse::
      :ref: jetstream.cli.subcommands.project.arg_parser
      :prog: jetstream project
      :nodefault:

project config
---------------

    ``jetstream project config``

    .. argparse::
       :ref: jetstream.cli.subcommands.project.config_arg_parser
       :prog: jetstream project config
       :nodefault:

project init
---------------

    ``jetstream project init``

    .. argparse::
       :ref: jetstream.cli.subcommands.project.init_arg_parser
       :prog: jetstream project init
       :nodefault:

project samples
---------------

    ``jetstream project samples``

    .. argparse::
       :ref: jetstream.cli.subcommands.project.samples_arg_parser
       :prog: jetstream project samples
       :nodefault:

project tasks
---------------

    ``jetstream project tasks``

    .. argparse::
       :ref: jetstream.cli.subcommands.project.tasks_arg_parser
       :prog: jetstream project tasks
       :nodefault:


project runs
---------------

    ``jetstream project runs``

    .. argparse::
       :ref: jetstream.cli.subcommands.project.runs_arg_parser
       :prog: jetstream project runs
       :nodefault:


legacy
-------

    ``jetstream legacy``

    .. argparse::
       :ref: jetstream.cli.subcommands.legacy.arg_parser
       :prog: jetstream legacy
       :nodefault:

report
-------

    ``jetstream report``

    .. argparse::
       :ref: jetstream.cli.subcommands.report.arg_parser
       :prog: jetstream report
       :nodefault:
