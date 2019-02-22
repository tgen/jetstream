"""Jetstream application settings

Settings are loaded in the module jetstream.settings with the help of
Confuse (https://github.com/sampsyo/confuse). The following is a description
of the loading process adapted from the Confuse docs:

Jetstream looks in a number of locations for application configuration files.
The locations are determined by the platform. For example, the first search
location on Unix-y systems is "$XDG_CONFIG_HOME/jetstream".

Here are the default search paths for each platform:

  OS X: ~/.config/jetstream and ~/Library/Application Support/jetstream
  Other Unix-y: $XDG_CONFIG_HOME/jetstream and ~/.config/jetstream
  Windows: %APPDATA%\jetstream where the APPDATA environment variable falls
    back to %HOME%\AppData\Roaming if undefined

Users can also add an override configuration directory with an environment
variable $JETSTREAMDIR.
"""
import json
import logging
import os
import jetstream

log = logging.getLogger(__name__)


def arg_parser(p):
    p.add_argument(
        '-f', '--full',
        action='store_true',
        help='Show the full settings values loaded from all sources'
    )


def main(args):
    log.debug(f'{__name__} {args}')
    full = jetstream.settings.flatten()

    if args.full:
        # hack to get rid of class info when serializing to yaml
        info = json.loads(json.dumps(full))
    else:
        info = {
            'Current setup': {
                'path': jetstream.settings.user_config_path(),
                'JETSTREAMDIR': os.environ.get('JETSTREAMDIR')
            }
        }

    print(__doc__)
    print(jetstream.utils.yaml_dumps(info))
