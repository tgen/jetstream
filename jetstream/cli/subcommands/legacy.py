import sys
import argparse
import logging

from jetstream import utils
from jetstream.legacy import config

log = logging.getLogger(__name__)


def build_parser():
    parser = argparse.ArgumentParser(
        prog='jetstream legacy',
        description="Convert legacy config files to YAML/JSON"
    )

    parser.add_argument('path')

    parser.add_argument('--format',
                        default='yaml',
                        choices=['yaml', 'json'])

    return parser


def main(args):
    parser = build_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    c = config.load(args.path)

    if args.format == 'yaml':
        utils.yaml.dump(c, stream=sys.stdout)
    elif args.format == 'json':
        print(utils.json.dumps(c, indent=4))
