import textwrap
import traceback
import logging
import jetstream


log = logging.getLogger(__name__)


class KvargsError(Exception):
    pass


def parse_key(k, separator=None):
    """Parses a kvarg key and returns Tuple(<loader>, <params>, <key>),
    where <loader> is the name used to fetch the loader function from the
    settings file, <params> are additional arguments passed to the loader fn,
    and <key> is the key that will be used to reference the loaded object. """
    if k.startswith('--'):
        k = k[2:]

    separator = separator or jetstream.settings['kvargs']['separator'].get(str)

    if not separator in k:
        err = f'Key "{k}" does not include separator "{separator}"'
        raise KvargsError(err)

    params, _, key = k.rpartition(separator)
    params = params.split(separator)
    loader = params.pop(0)

    return loader, params, key


def kvarg(k, v, separator=None):
    """Loads a kvarg pair and returns Tuple(key, obj).
    Kvargs follow the syntax where each item is delimited by separator:

        --<loader><params ... ><key> <value>

    The default separator is ':'. If the requested loader is not registered in
    the settings, or the loader function fails for the value, a KvargsError
    will be raised.
    """
    loader, params, key = parse_key(k, separator)

    try:
        lname = jetstream.settings['kvargs']['loaders'].get(dict)[loader]
    except KeyError:
        o = list(jetstream.settings['kvargs']['loaders'].all_contents())
        o = ', '.join(list(o))
        msg = f'Argument "{k} {v}": Unknown arg loader: "{loader}". Available' \
              f'arg loaders are: {o}.'
        raise KvargsError(msg) from None

    try:
        loader = jetstream.utils.dynamic_import(lname)
        obj = loader(v, *params)
    except Exception:
        tb = traceback.format_exc()
        msg = f'Error loading kvarg "{k} {v}":\n\n' \
              f'{textwrap.indent(tb, "  ")}\n' \
              f'Variable loading failed for  "{k} {v}", full traceback shown ' \
              f'above.'
        raise KvargsError(msg) from None

    return key, obj


def parse_kvargs(args, separator=None):
    """Parse arbitrary key-value arguments and return as a dictionary. """
    log.debug('Parsing kvargs: {}'.format(args))
    results = dict()
    iterator = iter(args)

    while 1:
        try:
            arg = next(iterator)
        except StopIteration:
            break

        if arg.startswith('--'):
            try:
                value = next(iterator)
            except StopIteration:
                msg = f'Argument {arg}: expected one following argument'
                raise KvargsError(msg) from None

            key, obj = kvarg(arg, value, separator=separator)
            results[key] = obj

    return results

