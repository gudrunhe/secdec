Installation
============

`pySecDec` should run fine with both, `python` 2.7 and `python` 3
on unix-like systems. It has been tested and developed on
`MacOS 10.11 (El Capitan)` and `openSUSE Tumbleweed`.

Before you install `pySecDec`, make sure that you have
recent versions of `numpy` (http://www.numpy.org/) and
`sympy` (http://www.sympy.org) installed.
Then download and upack the following tarball.

.. TODO: make `tarball` a download link of `complete_dist`

The tarball contains a distribution of `pySecDec` and
the additional dependencies listed :ref:`below <additional_cpp_dependencies>`.
Typing

.. code::

    $ make

should build all redistributed packages and display two commands
to be added to your ``.bashrc`` or ``.profile``.

.. note::
    Parallel build with ``make -j<number-of-cores>`` causes trouble
    on some systems. If ``make`` finished without a message starting
    with `Successfully built "pySecDec" and its dependencies`, try
    again without the ``-j`` option.

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
designed for `normaliz` versions ``3.0.0``, ``3.1.0``, and ``3.1.1``. We recommend
to set your ``$PATH`` such that the `normaliz` executable is found. Alternatively,
you can pass the path to the `normaliz` executable directly to the functions that
need it.

.. _additional_cpp_dependencies:

Additional Dependencies for Generated c++ Packages
--------------------------------------------------

The intended main usage of `pySecDec` is to make it write c++ packages using the functions
:func:`pySecDec.code_writer.make_package` and :func:`pySecDec.loop_integral.loop_package`.
In order to build these c++ packages, the following additional non-python-based libraries
and programs are required:

 * CUBA (http://www.feynarts.de/cuba/)
 * FORM (http://www.nikhef.nl/~form/)
 * gsl (http://www.gnu.org/software/gsl/)
 * SecDecUtil (part of `pySecDec`), depends on:

   * catch (https://github.com/philsquared/Catch)
   * fast-cpp-csv-parser (https://github.com/ben-strasser/fast-cpp-csv-parser)

These packages are redistributed with the `pySecDec` tarball; i.e. you don't have to install
any of them yourself.
