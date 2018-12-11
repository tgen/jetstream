#!/usr/bin/env python3
import os
import json
import shutil
import argparse
import logging
import subprocess

log = logging.getLogger(__name__)

try:
    with open(os.path.expanduser('~/.jetstream_dev')) as fp:
        settings = json.load(fp)
except Exception:
    log.exception('Unable to load dev settings json at ~/.jetstream_dev')
    settings = {}


def arg_parser(actions=None):
    parser = argparse.ArgumentParser(
        description='Deployment tool for Jetstream'
    )

    parser.add_argument('action', help='Action name', choices=actions)

    parser.add_argument('args', nargs=argparse.REMAINDER)

    return parser


def build_docs(args=None, open_index=True):
    subprocess.check_call(['make', 'html'], cwd='docs')

    if open_index:
        subprocess.call([
            'open',
            os.path.join('docs', '_build', 'html', 'index.html')
        ])


def publish_docs(args=None):
    build_docs(open_index=False)
    subprocess.check_call([
        'rsync',
        '--progress',
        '--delete',
        '-r',
        os.path.join('docs', '_build', 'html'),
        settings['deploy_docs_uri']
    ])

# This is not working. I think there is a bug in pip install when running
# inside of a python script, the binary path is not configured correctly
# def test(args=None):
#     try:
#         shutil.rmtree('venv')
#     except FileNotFoundError:
#         pass
#
#     subprocess.check_call('virtualenv venv', shell=True)
#     subprocess.check_call('source venv/bin/activate && which pip && pip install . && cat /Users/rrichholt/jetstream/venv/bin/jetstream', shell=True)
#     subprocess.check_call('source venv/bin/activate && which python && python3 -m unittest discover', shell=True)


def main(args=None):
    actions = {
        'build-docs': build_docs,
        'publish-docs': publish_docs,
    }

    parser = arg_parser(actions=actions)
    args = parser.parse_args(args)
    log.debug(args)

    actions[args.action](args.args)


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.DEBUG,
        format="[%(module)10s] %(asctime)s: %(message)s"
    )
    main()
