import logging


def verbose(self, message, *args, **kws):
    if self.isEnabledFor(5):
        self._log(5, message, args, **kws)


logging.addLevelName(5, 'VERBOSE')
logging.Logger.verbose = verbose

log = logging.getLogger('jetstream')
log.addHandler(logging.NullHandler())
log.setLevel(1)

basic_format = "{asctime}: {message}"
color_format = "[\033[4m\033[92m\U0001F335 {module:>10}\033[0m] " + basic_format
debug_format = "{levelname} {module}:{lineno} " + basic_format


def start_logging(*, format=basic_format, level=1):
    global log
    h = logging.StreamHandler()
    h.setFormatter(logging.Formatter(
        format,
        datefmt="%Y-%m-%d %H:%M:%S",
        style='{'
    ))
    h.setLevel(level)
    log.addHandler(h)
    return log
