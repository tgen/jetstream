# Command not found after install

Pip installs package scripts and executables automatically. In some causes you will need to update your [`$PATH`][linux_path] to include the directory where these scripts and executables are stored. 

[See this post][path_help] if you need more help


## MacOS

Changes should be added to `~/.bash_profile`. 

Installed system wide: the bin directory will be `/usr/local/bin/` and 

Installed with `--user`: the bin directory will be `~/Library/Python/<pythonversion>/bin`

## Linux

Changes should be added to `~/.bashrc`

Installed system wide: the bin directory will be `/usr/bin/`

Installed with `--user`: the bin directory will be `~/.local/bin/`



[path_help]: https://stackoverflow.com/questions/35898734/pip-installs-packages-successfully-but-executables-not-found-from-command-line
[linux_path]: http://www.linfo.org/path_env_var.html

