"""Shortcut to "jetstream run" with --build-only option enabled"""
from jetstream.cli.subcommands import run


def main(args=None):
    args.append('--build-only')
    return run.main(args)


if __name__ == '__main__':
    main()
