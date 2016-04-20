Installation
============

`pySecDec` should run fine with both, `python` 2.7 and `python` 3.
It has been tested and developed on `MacOS 10.11 (El Capitan)` and
`openSUSE Tumbleweed`. However, it should be platform independent
and also work on Windows.

.. _install_from_PyPI:

Installation from PyPI using `pip` (coming up soon)
----------------------------------------------------------------------

Installation is easiest using `pip` (https://pypi.python.org/pypi/pip).
`pip` automatically installs all dependencies
along with the desired package from the Python Package Index
(https://pypi.python.org/pypi).

.. code::

    $ pip install pySecDec

Manual Installation
-------------------

Before you manually install `pySecDec`, make sure that you have
recent versions of `numpy` (http://www.numpy.org/) and
`sympy` (http://www.sympy.org) installed.

To install `pySecDec`, open a shell in the top level directory (where
``setup.py`` is located) and type::

    $ python setup.py install

If you have `pip`, you should type

.. code::

    $ pip install .

instead, for the same reasons as mentioned :ref:`above <install_from_PyPI>`.

The Geomethod and Normaliz
--------------------------

.. note::
    If you are not urgently interested in using the
    :mod:`geometric decomposition <pySecDec.decomposition.geometric>`, you
    can ignore this section for the beginning. The instructions below are
    not essential for a `pySecDec` installation. You can still install
    `normaliz <https://www.normaliz.uni-osnabrueck.de/>`_
    **after** installing `pySecDec`. All but the
    :mod:`geometric decomposition <pySecDec.decomposition.geometric>`
    routines work without `normaliz <https://www.normaliz.uni-osnabrueck.de/>`_.

If you want to use the :mod:`geometric decomposition <pySecDec.decomposition.geometric>`
module, you need the
`normaliz <https://www.normaliz.uni-osnabrueck.de/>`_ [BIR]_ command line executable.
The :mod:`geometric decomposition <pySecDec.decomposition.geometric>` module is
designed for `normaliz <https://www.normaliz.uni-osnabrueck.de/>`_ version ``3.0.0``. We recommend to set your ``$PATH``
such that the `normaliz <https://www.normaliz.uni-osnabrueck.de/>`_ executable is found. Alternatively, you can pass the
path to the `normaliz <https://www.normaliz.uni-osnabrueck.de/>`_ executable directly to the functions that need it.

For Developers
--------------

`pip` offers an "editable" installation that can be triggered by::

    $ pip install -e /path/to/repository --user

This command causes `python` to load `pySecDec` directly from your local
copy of the repository. As a result, no reinstallation is required after
making changes in the source code.

`pySecDec` comes with a self test suite written in the `python unittest` framework.
The most convenient way to run all test is using `nose` (http://nose.readthedocs.org).
If `nose` is installed, just type::

    $ nosetests

in the source repository to run all tests. In order to check that the examples
given in the documentation are working, go to the ``doc`` subdirectory and type::

    $ make doctest

Also note the ``Makefile`` in the package's root directory that implements a
few common development tasks. You can list all available targets with the command
::

    $ make help
