pySecDec
========

[![Documentation Status](https://readthedocs.org/projects/secdec/badge/?version=latest)](http://secdec.readthedocs.io/en/latest/?badge=latest)

pySecDec is a toolbox for the calculation of dimensionally
regulated parameter integrals using the sector decomposition approach.

See [1811.11720](https://arxiv.org/abs/1811.11720)
and [1703.09692](https://arxiv.org/abs/1703.09692)
for description of the implementation;
[0803.4177](http://arxiv.org/abs/0803.4177)
and [hep-ph/0004013](http://arxiv.org/abs/hep-ph/0004013)
for description of sector decomposition; and
[1502.06595](http://arxiv.org/abs/1502.06595) for SecDec,
the predecessor of pySecDec.


Installation
============

pySecDec should work under Python version 3.6 or newer on the
usual Unix-like systems.

The latest release can be installed from [PyPI](https://pypi.org/project/pySecDec/)
by first (optionally) upgrading [pip](https://pypi.org/project/pip/):

    $ python3 -m pip install --user 'pip>=20.1'

and then running:

    $ python3 -m pip install --user --upgrade pySecDec

To install the latest development version on the other hand, run:

    $ git clone https://github.com/gudrunhe/secdec secdec
    $ cd secdec
    $ python3 -m pip install --user .

## The geometric method and Normaliz

If you want to use the geometric decomposition methods
(`decomposition_method='geometric'` or `'geometric_ku'`),
you need the `normaliz`  executable,
which can be downloaded from https://www.normaliz.uni-osnabrueck.de
[W. Bruns and B. Ichim and T. Römer and R. Sieg and C. Söger].
The geometric decomposition module is
designed for normaliz version 3 - currently versions
``3.3.0``, ``3.4.0``, ``3.5.4``, ``3.6.0``, ``3.6.2``, ``3.7.3``,
``3.7.4``, and ``3.8.1``
are known to work. We recommend to set your $PATH such that the
`normaliz` executable is found. Alternatively, you can pass the path to the `normaliz`
executable directly to the functions that need it.

## Additional dependencies for generated C++ packages

The intended main usage of pySecDec is to make it write C++ packages using the functions
`pySecDec.code_writer.make_package` and `pySecDec.loop_integral.loop_package`.
In order to build these C++ packages, the following additional non-python-based libraries
and programs are required:

 * CUBA (http://www.feynarts.de/cuba/)
 * FORM (http://www.nikhef.nl/~form/)
 * GSL (http://www.gnu.org/software/gsl/)
 * SecDecUtil (part of pySecDec), and its dependency:

   * Catch (https://github.com/philsquared/Catch)

The functions `pySecDec.code_writer.make_package` and
`pySecDec.loop_integral.loop_package` can use the external program
`dreadnaut` to find all sector symmetries and therefore reduce
the number of sectors:

 * NAUTY (http://pallini.di.uniroma1.it/)
[B. D. McKay and A. Piperno, Practical graph isomorphism, II, 2014, Journal of Symbolic Computation, 60, 94-112,
doi:10.1016/j.jsc.2013.09.003]

These packages are redistributed along with pySecDec itself,
and will be built automatically during pySecDec installation.


Basic Usage
===========

You can find the pySecDec manual over at https://secdec.readthedocs.io.
Additionally, the development repository contains a folder `examples` with examples of pySecDec usage.

A simple example of the evaluation of a loop integral with pySecDec is `box1L`.
This example computes a one-loop box with one off-shell leg (with off-shellness s1) and one internal massive line (with mass squared msq).

To run this example change into the `box1L` directory and run the commands:

    $ python3 generate_box1L.py
    $ make -C box1L
    $ python3 integrate_box1L.py

The file `generate_box1L.py` defines the loop integral and calls pySecDec to perform the sector decomposition.
Once run it produces the directory `box1L` which contains the code required to numerically evaluate the integral.
The `make` command builds this code and produces a library.
The file `integrate_box1L.py` loads the integral library and evaluates the integral for a specified numerical point.

This will print the result of the integral evaluated with Mandelstam invariants s=4.0, t=-0.75 and s1=1.25, msq=1.0:

    eps^-2: -0.142868356275422825 - 1.63596224151119965e-6*I +/- ( 0.00118022544307414272 + 0.000210769456586696187*I )
    eps^-1: 0.639405625715768089 + 1.34277036689902802e-6*I +/- ( 0.00650722394065588166 + 0.000971496627153705891*I )
    eps^0 : -0.425514350373418893 + 1.86892487760861536*I +/- ( 0.00706834403694714484 + 0.0186497890361357298*I )


## List of examples

 * `easy`: a simple parametric integral, described in the manual in Section 2.1.
 * `box1L`: a simple 1-loop, 4-point, 4-propagator integral, described in the manual Section 2.2.
 * `triangle2L`: a 2-loop, 3-point, 6-propagator diagram, also containing massive propagators.
 * `box2L_numerator`: a massless planar on-shell 2-loop, 4-point, 7-propagator box with a numerator, either defined as an inverse propagator `box2L_invprop.py` or in terms of contracted Lorentz vectors `box2L_contracted_tensor.py`.
 * `triangle3L`: a 2-loop, 3-point, 7-propagator integral, demonstrates that the symmetry finder can significantly reduce the number of sectors.
 * `elliptic2L_euclidean`: an integral known to contain elliptic functions, evaluated at a Euclidean phase-space point.
 * `elliptic2L_physical`: an integral known to contain elliptic functions, evaluated at a physical phase-space point.
 * `triangle2L_split`: a 2-loop, 3-point, 6-propagator integral without a Euclidean region due to special kinematics.
 * `hypergeo5F4`: a general dimensionally regulated parameter integral, corresponding to a Hypergeometric function 5F4.
 * `4photon1L_amplitude`: calculation of the 4-photon amplitude, showing how to use pySecDec as an integral library in a larger context.
 * `two_regulators`: an integral involving poles in two different regulators.
 * `userdefined_cpp`: a collection of examples demonstrating how to combine polynomials to be decomposed with other user-defined functions.


Development
===========

During development instead of full installation it is more
convenient to work directly from the checked out repository.
To make this work, first install the Python dependencies and
build the contributed software by running

    $ make dependencies
    $ make build

Then, define PYTHONPATH environment variable to point at the
root of the pySecDec sources.

The ``Makefile`` in the package's root directory also implements
other common development tasks.  You can list all available
targets with the command

    $ make help

`pySecDec` comes with a self test suite written in the `python unittest` framework.
The most convenient way to run all test is using `nose` (http://nose.readthedocs.org).
If `nose` is installed, type

    $ make check

in the source repository to run all tests. Developers should write test cases for
ALL functions they implement and make sure that ALL tests pass before uploading a
commit.

## Building the Documentation

To build the documentation of `pySecDec`, you need `sphinx` (http://www.sphinx-doc.org).
If `sphinx` is installed, the command

    $ make doc

generates the documentaion in `html` and in `pdf` format. Developers
should inline-document python functions and also keep the C++
part up to date. To generate HTML and PDF separately, use

    $ make doc-html

and

    $ make doc-pdf

Building the documentaion in pdf format requires an up-to-date installation of a latex
implementation. If you get an error about missing ".sty" file, do the following:

 1. If you are an administrator on your computer, try to install the missing latex packages
    with your favorite package manager. The MiKTeX or TeXLive implementations should contain
    all required packages.

 2. If you are not an administrator, first get the missing packages, e.g. from
    "http://www.ctan.org/". Collect the missing files in one or more directories
    and set the environment variable TEXINPUTS to these directories in the same
    way as the PATH variable is typically set.

## Making a PyPI release

To create a release on PyPI, first increment the version number
of pySecDec in `pyproject.toml`, and update `ChangeLog`. Then
create source distribution archive suitable for PyPI by running

    $ make dist

The archive will appear in `dist/pySecDec-<version>.tar.gz`. To
double-check it, list its contents:

    $ tar tzf dist/pySecDec-*.tar.gz

To upload this archive to PyPI follow the steps described in the
[Packaging Python Projects](https://packaging.python.org/tutorials/packaging-projects/) reference.
In short: once you've registered an account on pypi.org
and obtained a key, publishing a new release is done via
[Twine](https://pypi.org/project/twine/):

    $ twine check dist/pySecDec-<version>.tar.gz
    $ twine upload dist/pySecDec-<version>.tar.gz

### Uploading prebuilt releases

Note that along with source distribution (*.tar.gz) it is useful
to upload prebuilt version of it too, so that users would
not have to wait for the contributed software to build. To
this end, one should upload `dist/pySecDec-*.whl` files
too. The caveat is that these must be built inside one of the
[manylinux](https://github.com/pypa/manylinux) Docker images:
this guarantees that the prebuilt libraries and programs will
work on a wide range of Linux systems. Once this is done, double
check that the built distribution has the correct platform tag
by using [auditwheel](https://pypi.org/project/auditwheel/):

    $ auditwheel show dist/pySecDec-<version>-<tag>.whl
    
And then upload:

    $ twine check dist/pySecDec-<version>-<tag>.whl
    $ twine upload dist/pySecDec-<version>-<tag>.whl
