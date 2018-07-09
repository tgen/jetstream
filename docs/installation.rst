Installation
============

Python3
-------

.. note::
    **TGen users on Dback can load the latest version of Python with**
    ``module load python`` **and skip to the next step.**


Installation requires Python3 and has only been tested on MacOS and
Linux. Python installation guides are available from the `Hitchiker’s
Guide to Python`_.



Install with Pip
----------------

`Pip`_ is the official package manager for Python. It simplifies
installation of Python packages by automatically downloading and
installing dependencies. Choose one of the commands below to install the
latest version of Jetstream into your user package library:

Install via HTTPS (requires entering Github username and password)

.. code:: sh

   pip3 install --upgrade --user git+https://github.com/tgen/jetstream.git@master

Or, install via SSH (requires SSH keys to be configured in Github
profile)

.. code:: sh

   pip3 install --upgrade --user git+ssh://git@github.com/tgen/jetstream.git@master

Verify Install
--------------

The command-line utilites should now be available. Try them out:

.. code:: sh

   jetstream --help

If that prints a help message, you’re ready to go! If you get an error
instead, it’s likely that the Python bin directory is not added to your
path. `See troubleshooting`_ for more help.

.. _Hitchiker’s Guide to Python: http://docs.python-guide.org/en/latest/starting/installation/
.. _Pip: https://pip.pypa.io
.. _See troubleshooting: troubleshooting.rst
