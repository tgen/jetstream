Troubleshooting
===============

Command not found after install
-------------------------------

Pip installs package scripts and executables automatically. In some
causes you will need to update your `PATH`_ to include the
directory where these scripts and executables are stored.

`See this post`_ if you need more help

- **MacOS**:

    Changes should be added to ``~/.bash_profile``.

    Installed system-wide: the bin directory will be ``/usr/local/bin/``

    Installed with ``--user``: the bin directory will be
    ``~/Library/Python/<pythonversion>/bin``

- **Linux**:

    Changes should be added to ``~/.bashrc``

    Installed system wide: the bin directory will be ``/usr/bin/``

    Installed with ``--user``: the bin directory will be ``~/.local/bin/``

.. _PATH: http://www.linfo.org/path_env_var.html
.. _See this post: https://stackoverflow.com/questions/35898734/pip-installs-packages-successfully-but-executables-not-found-from-command-line
