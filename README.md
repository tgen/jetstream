# Installation

## Install Python3 and Pip

> TGen users on Dback can load the latest version of Python with `module load python` and skip to the next step.

This is a Python package and requires Python3. Installation guides for Mac/Windows/Linux are available from the [Hitchiker's Guide to Python][install_help] After Python3 is installed, you can install Jetstream with Pip (next step).


## Install Jetstream with Pip

Install via HTTPS (requires entering Github username and password)

```shell
pip3 install --upgrade --user git+https://github.com/tgen/jetstream.git@master
```

Install via SSH (requires SSH keys to be configured in Github profile)

```shell
pip3 install --upgrade --user git+ssh://git@github.com/tgen/jetstream.git@master
```

# Docs

Docs are being hosted here http://dback-login3.tgen.org:8082 , they're still largely a work-in-progress:



[install_help]: http://docs.python-guide.org/en/latest/starting/installation/
[path_help]: https://stackoverflow.com/questions/35898734/pip-installs-packages-successfully-but-executables-not-found-from-command-line
