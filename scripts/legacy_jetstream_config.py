import sys
import argparse
import jetstream
from jetstream.legacy.config import load, explode

log = jetstream.log


def arg_parser():
    parser = argparse.ArgumentParser(
        prog='jetstream legacy',
        description="Convert legacy config files to YAML/JSON"
    )

    parser.add_argument('path', help='Path to a legacy .config file.')

    parser.add_argument('--json', dest='format', action='store_const',
                        const='json', help='Output JSON')

    parser.add_argument('--yaml', dest='format', action='store_const',
                        const='yaml', help='Output YAML')

    parser.add_argument('--explode', dest='format', action='store_const',
                        const='explode', help='Explode into multiple files')

    parser.add_argument('--format',
                        default='yaml',
                        choices=['yaml', 'json', 'explode'])

    return parser



def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    c = load(args.path)

    if args.format == 'yaml':
        jetstream.utils.yaml_dump(c, stream=sys.stdout)

    elif args.format == 'json':
        jetstream.utils.json_dump(c, fp=sys.stdout)

    elif args.format == 'explode':
        explode(c)


if __name__ == '__main__':
    main()
