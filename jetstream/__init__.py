"""Jetstream is a collection of tools for automating workflows at TGen."""
import os

__version__ = '0.1.0a1'
__author__ = 'Ryan Richholt'
# TODO finish these attributes

# This prevents numpy from starting a million threads when imported. The
# graph library, networkx, uses scipy. TODO switch to another graph lib
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
