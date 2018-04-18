import os
import sys
import argparse
import subprocess
import logging

log = logging.getLogger(__name__)

def arg_parser():
    parser = argparse.ArgumentParser(
        description="Helper for launching jetstream pipeline components. All "
                    "arguments following first positional <cmd> will be passed"
                    "to command. Any arguments meant for this utility must come"
                    "before first positional"
    )

    parser.add_argument('cmd', nargs=argparse.REMAINDER)

    parser.add_argument('-r', '--release')
    # TODO how would a default release work?

    return parser


def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    release_dir = os.path.join(os.environ['JETSTREAM_RELEASES_DIR'], args.release)
    if not os.path.exists(release_dir):
        raise ValueError('Release not found: {}'.format(release_dir))

    os.environ["PATH"] += os.pathsep + release_dir

    log.debug(os.environ["PATH"])

    log.debug(args.cmd)
    p = subprocess.Popen(args.cmd)
    sys.exit(p.wait())
