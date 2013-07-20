"""
Configuration utilities
"""
import os
from configobj import ConfigObj

# Modified from https://github.com/tonysyu/mpltools/blob/master/setup.py,
# accessed 6/30/2013

__all__ = ['iter_paths', 'read', 'config']


def iter_paths(config_paths):
    for path in config_paths:
        path = os.path.expanduser(path)

        if not os.path.exists(path):
            continue

        yield read(path)


def read(path):
    """Return dict-like object of config parameters from file path.
    """
    return ConfigObj(path)


# Set fcres specific properties (i.e., not matplotlib properties).
config = {}
pkgdir = os.path.abspath(os.path.dirname(__file__))
for cfg in iter_paths([os.path.join(pkgdir, 'fcres.ini'),
                       '~/.fcres.ini',
                       './fcres.ini']):
    config.update(cfg)
