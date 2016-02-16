Installation
============

`pySecDec` should run fine with both, `python` 2.7 and `python` 3.
It has been tested and developed on `MacOS 10.11 (El Capitan)` and
`openSUSE 13.2 (Harlequin)`. However, it should be platform independent
and also work in Windows.

Installation from PyPI using `pip` (recommended, but not possible yet)
----------------------------------------------------------------------

Installation is easiest using `pip` (https://pypi.python.org/pypi/pip).
`pip` automatically installs all dependencies
along with the desired package from the Python Package Index
(https://pypi.python.org/pypi).
::

    $ pip install pySecDec

Manual Installation
-------------------

Before you manually install `pySecDec`, make sure that you have
recent versions of `numpy` (http://www.numpy.org/) and
`sympy` (http://www.sympy.org) installed.

To install `pySecDec`, open a shell in the source repository and type::

    $ python setup.py install

For Developers
--------------

`pip` offers an "editable" installation, that can be triggered by::

    $ pip install -e /path/to/repository --user

This command causes `python` to load `pySecDec` directly from your local
copy of the repository. As a result, no reinstallation is required after
making changes in the source code.
