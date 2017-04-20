"""
Code Writer
-----------

This module collects routines to create a c++ library.

Make Package
~~~~~~~~~~~~

This is the main function of `pySecDec`.

.. autofunction:: pySecDec.code_writer.make_package

Template Parser
~~~~~~~~~~~~~~~

.. automodule:: pySecDec.code_writer.template_parser
    :members:

.. include:: generated_cpp.txt


"""

from . import template_parser
from .make_package import make_package
