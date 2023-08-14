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

`pySecDec` works under Python version 3.8 or newer on
unix-like systems.  The latest release can be installed from
`PyPI`_ by first (optionally) upgrading `pip`_:

.. code::

    $ python3 -m pip install --user 'pip>=20.1'

and then running:

.. code::

   $ python3 -m pip install --user --upgrade pySecDec

This command will install the prebuild version of `pySecDec` if it
is available; if not, then the dependencies will be compiled from
source (this might take a while). One can also force building
from source like this:

.. code::

   $ python3 -m pip install --user --upgrade --no-binary :all: pySecDec

.. _PyPI: https://pypi.org/project/pySecDec/
.. _pip: https://pypi.org/project/pip/

.. _installation_neato:

Drawing Feynman Diagrams with `neato`
-------------------------------------

The :func:`~pySecDec.loop_integral.draw.plot_diagram` function draws a Feynman diagram using the command line tool `neato`. 
It is automatically called when generating an integral library using the :func:`~pySecDec.loop_integral.loop_package` function with an integral of type  :class:`~pySecDec.loop_integral.LoopIntegralFromGraph` and will issue a warning if `neato` is not available. 

The warning can be safely ignored if you are not interested in the drawing.
Alternatively, you must manually install `neato` which is part of the `graphviz` package.
It is available in many package repositories and at http://www.graphviz.org.

.. _additional_cpp_dependencies:

Additional Dependencies
-----------------------

`pySecDec` and the integration libraries it produces depend
on multiple third-party non-Python packages, all of which are
contained in `pySecDecContrib` and will be automatally built
during the normal installation procedure. These packages are:

 * QMC (https://github.com/mppmu/qmc), used for the :mod:`Qmc<pySecDec.integral_interface.Qmc>` integrator.
 * CUBA (http://www.feynarts.de/cuba/), used for :mod:`Vegas<pySecDec.integral_interface.Vegas>`, :mod:`Suave<pySecDec.integral_interface.Suave>`, :mod:`Divonne<pySecDec.integral_interface.Divonne>`, and :mod:`Cuhre<pySecDec.integral_interface.Cuhre>` integrators.
 * GSL (http://www.gnu.org/software/gsl/), used for the :mod:`CQuad<pySecDec.integral_interface.CQuad>` integrator.
 * FORM (http://www.nikhef.nl/~form/), used to optimize the integrands.
 * Nauty and Traces (http://pallini.di.uniroma1.it/), used by :func:`pySecDec.make_package` to find symmetries between sectors (if `use_dreadnaut` is set to `True`).
 * Normaliz (https://www.normaliz.uni-osnabrueck.de), used by the :mod:`geometric decomposition <pySecDec.decomposition.geometric>` module.
 * Catch (https://github.com/philsquared/Catch) used by :ref:`SedDecUtil<chapter_secdecutil>` for unit testing.
