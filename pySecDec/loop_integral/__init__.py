"""
Loop Integral
-------------

This module defines routines to Feynman parametrize a
loop integral and build a c++ package that numerically
integrates over the sector decomposed integrand.

Feynman Parametrization
~~~~~~~~~~~~~~~~~~~~~~~

Routines to Feynman parametrize a loop integral.

.. autoclass:: pySecDec.loop_integral.LoopIntegral
.. autoclass:: pySecDec.loop_integral.LoopIntegralFromGraph
.. autoclass:: pySecDec.loop_integral.LoopIntegralFromPropagators

Loop Package
~~~~~~~~~~~~

This module contains the function that generates a c++ package.

.. autofunction:: pySecDec.loop_integral.loop_package

Drawing Feynman Diagrams
~~~~~~~~~~~~~~~~~~~~~~~~

Use the following function to draw Feynman diagrams.

.. autofunction:: pySecDec.loop_integral.draw.plot_diagram

Loop Regions
~~~~~~~~~~~~

Applies the expansion by regions method to a loop integral.

.. autofunction:: pySecDec.loop_integral.loop_regions

"""

from .common import LoopIntegral
from .from_graph import LoopIntegralFromGraph
from .from_propagators import LoopIntegralFromPropagators
from .loop_package import loop_package, LoopPackage
from .draw import plot_diagram
from .loop_regions import loop_regions
