from . import start, resume, config, launch, plugins
__all__ = ["start", "resume", "config", "launch", "plugins"]


# TODO add subparsers for plugins module
# Need to be able to sync plugins library, update, maybe add/remove


# TODO add subparsers for other use cases of the app
# For example, pure shell script plugins can use other
# cli commands to make it easier for accessing project
# data or transforming genomic data
