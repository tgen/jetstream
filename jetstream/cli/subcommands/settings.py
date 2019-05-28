"""Jetstream Application Settings

Example:
jetstream settings --create --backend slurm --pipelines-home ~/pipelines

Use --full to see all current settings, or --create to initialize a new
config file. Use settings -h/--help to for more info on creating settings
files.
"""
import json
import logging
import os
import jetstream

log = logging.getLogger(__name__)
template = """# Jetstream Common User Settings
backend: {backend}
pipelines:
  home: {pipelines_home}

"""

# Details:
#
# Settings are loaded in the module jetstream.settings with the help of
# Confuse (https://github.com/sampsyo/confuse). The following is a description
# of the loading process adapted from the Confuse docs:
#
# Jetstream looks in a number of locations for application configuration files.
# The locations are determined by the platform. For example, the first search
# location on Unix-y systems is "$XDG_CONFIG_HOME/jetstream". The environment
# variable $JETSTREAMDIR can override this search path.
#
# Here are the default search paths for each platform:
#
#   OS X: ~/.config/jetstream and ~/Library/Application Support/jetstream
#   Other Unix-y: $XDG_CONFIG_HOME/jetstream and ~/.config/jetstream
#   Windows: %APPDATA%\jetstream where the APPDATA environment variable falls
#     back to %HOME%\AppData\Roaming if undefined


def arg_parser(parser):
    parser.add_argument(
        '--full',
        action='store_true',
        help='Show the full settings values loaded from all sources'
    )

    create = parser.add_argument_group('Create a new settings file')

    create.add_argument(
        '--create',
        action='store_true',
        help='Initialize an example config file at the user config path'
    )

    create.add_argument(
        '--backend',
        default='local',
        help='Backend to use when initializing a user settings file'
    )

    create.add_argument(
        '--pipelines-home',
        default='null',
        help='Pipelines home directory to use when initializing a user'
             'settings file.'
    )

    create.add_argument(
        '--overwrite-existing',
        action='store_true',
        help='Ignore FileExistsError when creating a settings file'
    )


def main(args):
    log.debug(f'{__name__} {args}')
    path = jetstream.settings.user_config_path()

    if args.full:
        full = jetstream.settings.flatten()
        print(jetstream.utils.yaml_dumps(json.loads(json.dumps(full))))
        return

    if args.create:
        if os.path.exists(path) and not args.overwrite_existing:
            raise FileExistsError(path)
        else:
            with open(path, 'w') as fp:
                settings = template.format(
                    backend=args.backend,
                    pipelines_home=args.pipelines_home
                )
                fp.write(settings)
            log.info(f'Created settings file at: {path}')
    else:
        info = f"Current setup:\npath: {path}\nexists: {os.path.exists(path)}"
        print(__doc__)
        print(info)
