import argparse
import logging
from jetstream import utils
from jetstream.core.legacy import config

log = logging.getLogger(__name__)


def build_parser():
    parser = argparse.ArgumentParser(
        description="Convert legacy config files to YAML"
    )

    parser.add_argument('path')

    return parser


def main(args):
    parser = build_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    c = config.load(args.path)
    print(utils.yaml_dumps(c))
