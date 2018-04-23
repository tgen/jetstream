#!/usr/bin/env python3
"""List data records from table"""
import json
import argparse
import jetstream


def arg_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument('key', help='Which config data to list records from')

    return parser


def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)

    p = jetstream.Project()

    for record in p.config[args.key]:
        print(json.dumps(record))


if __name__ == "__main__":
    main()
