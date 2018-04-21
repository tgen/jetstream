import argparse
import subprocess
import sys

import pkg_resources

import jetstream


def arg_parser():
    parser = argparse.ArgumentParser(
        description='Run a jetstream pipeline'
    )

    parser.add_argument('workflow')

    parser.add_argument('-v', '--version', action='version',
                        version=jetstream.__version__)

    return parser


def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)

    resource = pkg_resources.resource_filename(
        'jetstream_pipelines',
        'templates/{}.yaml'.format(args.workflow)
    )

    cmd = [
        'jetstream',
        'run',
        resource
    ]

    sys.exit(subprocess.call(cmd))



if __name__ == '__main__':
    main()
