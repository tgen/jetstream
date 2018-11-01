Command-Line
============


    ``jetstream``

    All command-line features can be accessed with ``jetstream <command>``
    Arguments described in this section can be used with any command.

    .. argparse::
        :ref: jetstream.cli.jetstream.arg_parser
        :prog: jetstream


init
------

    ``jetstream init``

    .. argparse::
        :ref: jetstream.cli.subcommands.init.arg_parser
        :prog: jetstream init
        :nodefault:


run
----------

    ``jetstream run``


    .. argparse::
        :ref: jetstream.cli.subcommands.run.arg_parser
        :prog: jetstream run
        :nodefault:


project
--------

    ``jetstream project``


    .. argparse::
      :ref: jetstream.cli.subcommands.project.arg_parser
      :prog: jetstream project
      :nodefault:


draw
------

    ``jetstream draw``

    .. argparse::
        :ref: jetstream.cli.subcommands.draw.arg_parser
        :prog: jetstream draw
        :nodefault:


build
-------

    ``jetstream build``

    .. argparse::
        :ref: jetstream.cli.subcommands.build.arg_parser
        :prog: jetstream build
        :nodefault:
