"""Tgen Rsync will recognize qualified URIs and initiate remote transfers

This will require that ssh keys are setup for the machines listed in
settings['datamovers']

"""
import argparse
import logging
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import urlparse, urlunparse
import jetstream
from jetstream.utils import data_transfer

log = logging.getLogger(__name__)


def resolve_location(value):
    """Given a file path, convert to valid uri"""
    if isinstance(value, str):
        uri = urlparse(value, scheme='file')
    else:
        uri = value

    try:
        wrangler = data_transfer.RsyncLimitedWrangler(
            hosts=jetstream.settings['data_movers'],
            max_rsyncs=jetstream.settings.get('max_rsyncs', 3)
        )
    except KeyError:
        log.warning('No data movers have been configured in settings profile.')
        return uri

    if uri.scheme == 'file':
        if uri.netloc == 'isilon.tgen.org':
            new_host = wrangler.next()
            uri = uri._replace(netloc=new_host)

    return uri


def rsync(src, dest, *args):
    src_uri = resolve_location(src)
    dest_uri = resolve_location(dest)

    log.info('Source: {}'.format(urlunparse((src_uri))))
    log.info('Destination: {}'.format(urlunparse(dest_uri)))

    ss = src_uri.scheme
    ds = dest_uri.scheme

    if ss == 'file' and ds == 'file':
        return data_transfer.rsync(src_uri, dest_uri, *args)
    else:
        raise NotImplementedError('No method for: {} -> {}'.format(ss, ds))


def main(args=None):
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-workers', type=int, default=1)
    parser.add_argument('paths', nargs='+')
    args, unknown = parser.parse_known_args(args)

    if len(args.paths) < 2:
       raise ValueError('No destination given')
    else:
        src_paths = args.paths[:-1]
        dest_path = args.paths[-1]

    executor = ThreadPoolExecutor(max_workers=args.max_workers)

    with executor:
        for src in src_paths:
            executor.submit(rsync, src, dest_path, *unknown)


