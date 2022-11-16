Installation
============

Obtain a Compiler
-----------------

`pySecDec` works on Unix-like systems, specifically Linux and macOS. It requires a working c++ compiler and a Python 3 (>=3.8) installation.

On Linux systems, users can install compilers (we recommend the latest GCC or Clang compiler) and Python 3 using the package manager provided with their Linux distribution (usually one of : `apk`, `apt-get`, `apt`, `yum`).

On macOS systems, users can install a compiler via the App Store:

* Download `Xcode` via the App Store
* Open a terminal window and enter :code:`xcode-select --install` then follow the prompts to install command line tools.

This procedure will make the `clang` compiler available on the system. Python 3 is shipped with macOS, you can check the version available on your system using :code:`python3 --version`. Later versions of Python can be installed on macOS using third-party package managers (e.g. the Homebrew package manager).


Download the Program and Install
--------------------------------

`pySecDec` works under Python version 3.6 or newer on
unix-like systems.  The latest release can be installed from
`PyPI`_ by first (optionally) upgrading `pip`_:

.. code::

    $ python3 -m pip install --user 'pip>=20.1'

and then running:

.. code::

   $ python3 -m pip install --user --upgrade pySecDec

.. _PyPI: https://pypi.org/project/pySecDec/
.. _pip: https://pypi.org/project/pip/

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

.. note::
    The following packages are redistributed with the `pySecDec` tarball; i.e. you don't have 
    to install any of them yourself.

The intended main usage of `pySecDec` is to make it write c++ packages using the functions
:func:`pySecDec.code_writer.make_package` and :func:`pySecDec.loop_integral.loop_package`.
In order to build these c++ packages, the following additional non-python-based libraries
and programs are used:

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

The :mod:`geometric decomposition <pySecDec.decomposition.geometric>`
module depends on the `normaliz` [BIR]_ command line executable:

 * Normaliz (https://www.normaliz.uni-osnabrueck.de)

These packages are redistributed along with pySecDec itself,
and will be built automatically during pySecDec installation.
