def _import_formats():
    from os import listdir
    from os.path import dirname, join, isfile, splitext, basename

    __all__ = []
    dn = dirname(__file__)
    for f in listdir(dn):
        fullpath = join(dn, f)
        if (isfile(fullpath) and f.endswith(('.py', '.pyc')) and not
        f.startswith('_')):
            name = splitext(basename(f))[0]
            __all__.append(name)
    return __all__

__all__ = _import_formats()
del _import_formats
