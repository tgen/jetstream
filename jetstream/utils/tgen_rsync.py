"""Tgen Rsync will recognize qualified URIs and initiate remote transfers

This will require that ssh keys are setup for the machines listed in
settings['datamovers']

"""
import argparse
import logging
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
            max_rsyncs=jetstream.settings.get('max_rsyncs', 4)
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
    parser.add_argument('src')
    parser.add_argument('dest')
    parser.add_argument('extra_args', nargs=argparse.REMAINDER,
                        help='Additional arguments will be passed to rsync')
    args = parser.parse_args(args)


    return rsync(args.src, args.dest, *args.extra_args)

