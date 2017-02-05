'''
Decomposition
-------------

The core of sector decomposition. This module implements
the actual decomposition routines.

Common
~~~~~~

This module collects routines that are used by
multiple decompition modules.

.. autoclass:: pySecDec.decomposition.Sector
.. autofunction:: pySecDec.decomposition.squash_symmetry_redundant_sectors_sort
.. autofunction:: pySecDec.decomposition.squash_symmetry_redundant_sectors_dreadnaut

Iterative
~~~~~~~~~

.. automodule:: pySecDec.decomposition.iterative
    :members:

Geometric
~~~~~~~~~

.. automodule:: pySecDec.decomposition.geometric
    :members:

Splitting
~~~~~~~~~

.. automodule:: pySecDec.decomposition.splitting
    :members:

'''

from . import iterative, geometric, splitting
from .common import *
