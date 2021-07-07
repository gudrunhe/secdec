"""
Code Writer
-----------

This module collects routines to create a c++ library.

Make Package
~~~~~~~~~~~~

This is the main function of `pySecDec`.

.. autofunction:: pySecDec.code_writer.make_package

Sum Package
~~~~~~~~~~~

Computing weighted sums of integrals, e.g. amplitudes.

.. automodule:: pySecDec.code_writer.sum_package
    :members:

Template Parser
~~~~~~~~~~~~~~~

.. automodule:: pySecDec.code_writer.template_parser
    :members:

.. include:: generated_cpp.txt

"""

from . import template_parser
from .make_package import make_package, MakePackage
from .sum_package import sum_package, Coefficient
