Installation
============

Download the Program and Install
--------------------------------

`pySecDec` should run fine with both, `python` 2.7 and `python` 3
on unix-like systems.

Before you install `pySecDec`, make sure that you have
recent versions of `numpy` (http://www.numpy.org/) and
`sympy` (http://www.sympy.org/) installed. The version of `sympy`
should be 0.7.6 or higher, the version of `numpy` should be 1.6 or higher.
Type

.. code::

   $ python -c "import numpy"
   $ python -c "import sympy"

to check for their availability.

In case either `numpy` or `sympy` are missing on your machine,
it is easiest to install them from your package repository. Alternatively,
and in particular if you do not have administrator rights,
`pip` (https://pip.pypa.io/en/stable/) may be used to perform
the installation.

To install `pySecDec` download and upack the tarball from http://secdec.hepforge.org/.
The tarball contains a distribution of `pySecDec` and the additional dependencies
listed :ref:`below <additional_cpp_dependencies>`. Typing

.. code::

    $ make

should build all redistributed packages and display two commands
to be added to your ``.bashrc`` or ``.profile``.

.. _installation_normaliz:

The Geomethod and Normaliz
--------------------------

.. note::
    If you are not urgently interested in using the
    :mod:`geometric decomposition <pySecDec.decomposition.geometric>`, you
    can ignore this section for the beginning. The instructions below are
    not essential for a `pySecDec` installation. You can still install
    `normaliz` **after** installing `pySecDec`. All but the
    :mod:`geometric decomposition <pySecDec.decomposition.geometric>`
    routines work without `normaliz`.

If you want to use the :mod:`geometric decomposition <pySecDec.decomposition.geometric>`
module, you need the `normaliz` [BIR]_ command line executable.
The :mod:`geometric decomposition <pySecDec.decomposition.geometric>` module is
designed for `normaliz` version 3 - currently versions ``3.0.0``, ``3.1.0``, ``3.1.1``,
``3.3.0``, ``3.4.0``, and ``3.5.4`` are known to work. We recommend to set your ``$PATH`` such that the
`normaliz` executable is found. Alternatively, you can pass the path to the `normaliz`
executable directly to the functions that need it.

.. _installation_neato:

Drawing Feynman Diagrams with `neato`
-------------------------------------

In order to use :func:`~pySecDec.loop_integral.draw.plot_diagram`, the command line tool
`neato` must be available. The function :func:`~pySecDec.loop_integral.loop_package` tries
to call :func:`~pySecDec.loop_integral.draw.plot_diagram` if given a
:class:`~pySecDec.loop_integral.LoopIntegralFromGraph` and issues a warning on failure. That
warning can be safely ignored if you are not interested in the drawing.

`neato` is part of the `graphviz` package. It is available in many package repositories and at
http://www.graphviz.org.

.. _additional_cpp_dependencies:

Additional Dependencies for Generated c++ Packages
--------------------------------------------------

The intended main usage of `pySecDec` is to make it write c++ packages using the functions
:func:`pySecDec.code_writer.make_package` and :func:`pySecDec.loop_integral.loop_package`.
In order to build these c++ packages, the following additional non-python-based libraries
and programs are required:

 * CUBA (http://www.feynarts.de/cuba/)
 * QMC (https://github.com/mppmu/qmc)
 * FORM (http://www.nikhef.nl/~form/)
 * SecDecUtil (part of `pySecDec`, see :ref:`SedDecUtil<chapter_secdecutil>`), depends on:

   * catch (https://github.com/philsquared/Catch)
   * gsl (http://www.gnu.org/software/gsl/)

The functions :func:`pySecDec.code_writer.make_package` and :func:`pySecDec.loop_integral.loop_package`
can use the external program `nauty` [MP+14]_ to find all sector symmetries and therefore reduce the number of
sectors:

 * NAUTY (http://pallini.di.uniroma1.it/)

These packages are redistributed with the `pySecDec` tarball; i.e. you don't have to install
any of them yourself.
