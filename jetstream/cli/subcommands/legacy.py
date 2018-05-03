import os
import sys
import csv
import argparse
import logging
import jetstream

log = logging.getLogger(__name__)


def build_parser():
    parser = argparse.ArgumentParser(
        prog='jetstream legacy',
        description="Convert legacy config files to YAML/JSON"
    )

    parser.add_argument('path')

    parser.add_argument('--json', dest='format', action='store_const',
                        const='json')

    parser.add_argument('--yaml', dest='format', action='store_const',
                        const='yaml')

    parser.add_argument('--explode', dest='format', action='store_const',
                        const='explode')

    parser.add_argument('--format',
                        default='yaml',
                        choices=['yaml', 'json', 'explode'])

    return parser


def records_to_csv(records, outpath):
    log.debug('to csv: {}\n{}'.format(outpath, records))

    keys = set()
    for row in records:
        keys = keys.union(set(row.keys()))

    log.debug('Found keys: {}'.format(keys))

    with open(outpath, 'w') as fp:
        dw = csv.DictWriter(fp, keys)
        dw.writeheader()
        dw.writerows(records)


def main(args):
    parser = build_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    c = jetstream.legacy.config.load(args.path)

    if args.format == 'yaml':
        jetstream.utils.yaml.dump(c, stream=sys.stdout)

    elif args.format == 'json':
        jetstream.utils.json.dump(c, fp=sys.stdout)

    elif args.format == 'explode':
        outdir = os.path.splitext(os.path.basename(args.path))[0]
        os.mkdir(outdir)

        log.critical('Exploding into: {}'.format(outdir))
        records_to_csv(c['samples'], os.path.join(outdir, 'samples.csv'))
        records_to_csv(c['data'], os.path.join(outdir, 'data.csv'))

        with open(os.path.join(outdir, 'run_parameters.yaml'), 'w') as fp:
            jetstream.utils.yaml.dump(c['run_parameters'], stream=fp)
