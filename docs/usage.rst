Command-Line
============


    ``jetstream``

    All command-line features can be accessed with ``jetstream <command>``
    Arguments described in this section can be used with any command.

    .. argparse::
        :ref: jetstream.cli.jetstream.arg_parser
        :prog: jetstream

run
----------

    ``jetstream run``


    .. argparse::
        :ref: jetstream.cli.subcommands.run.arg_parser
        :prog: jetstream run
        :nodefault:

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


project tasks
---------------

    ``jetstream project tasks``

    .. argparse::
       :ref: jetstream.cli.subcommands.project.tasks_arg_parser
       :prog: jetstream project tasks
       :nodefault:


project history
---------------

    ``jetstream project history``

    .. argparse::
       :ref: jetstream.cli.subcommands.project.history_arg_parser
       :prog: jetstream project history
       :nodefault:

