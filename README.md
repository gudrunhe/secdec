pySecDec
========

[![Documentation Status](https://readthedocs.org/projects/secdec/badge/?version=latest)](http://secdec.readthedocs.io/en/latest/?badge=latest)

pySecDec is a toolbox for the calculation of dimensionally
regulated parameter integrals using the sector decomposition approach.

See [2305.19768], [2108.10807], [1811.11720], and [1703.09692] for description of the
implementation; [0803.4177] and [hep-ph/0004013] for description
of sector decomposition; and [1502.06595] for SecDec, the
predecessor of pySecDec.

[2305.19768]: https://arxiv.org/abs/2305.19768
[2108.10807]: https://arxiv.org/abs/2108.10807
[1811.11720]: https://arxiv.org/abs/1811.11720
[1811.11720]: https://arxiv.org/abs/1811.11720
[1703.09692]: https://arxiv.org/abs/1703.09692
[0803.4177]: http://arxiv.org/abs/0803.4177
[hep-ph/0004013]: http://arxiv.org/abs/hep-ph/0004013
[1502.06595]: http://arxiv.org/abs/1502.06595

Installation
============

pySecDec should work under Python version 3.8 or newer on the
usual Unix-like systems.

The latest release can be installed from [PyPI] by first
upgrading [pip]:

    $ python3 -m pip install --user 'pip>=20.1'

and then running:

    $ python3 -m pip install --user --upgrade pySecDec

This command will install the prebuild version of `pySecDec` if it
is available; if not, then the dependencies will be compiled from
source (this might take a while). One can also force building
from source like this:

    $ python3 -m pip install --user --upgrade --no-binary :all: pySecDec

[pypi]: https://pypi.org/project/pySecDec/
[pip]: https://pypi.org/project/pip/

## Included dependencies

`pySecDec` and the integration libraries it produces depend
on multiple third-party non-Python packages, all of which are
contained in `pySecDecContrib/` folder and will be automatally
built during the normal installation procedure. These packages
are:

 * QMC (https://github.com/mppmu/qmc), used for the `Qmc` integrator.
 * CUBA (http://www.feynarts.de/cuba/), used for `Vegas`, `Suave`, `Divonne`, and `Cuhre` integrators.
 * GSL (http://www.gnu.org/software/gsl/), used for the `CQuad` integrator.
 * FORM (http://www.nikhef.nl/~form/), used to optimize the integrands.
 * zlib (http://zlib.net/), used by FORM.
 * Nauty and Traces (http://pallini.di.uniroma1.it/), used by `pySecDec.make_package` to find symmetries between sectors (if `use_dreadnaut` is set to `True`).
 * Normaliz (https://www.normaliz.uni-osnabrueck.de), used by the geometric decomposition module (`decomposition_method` set to `'geometric'` or `'geometric_ku'`).
 * Catch (https://github.com/philsquared/Catch) used by SedDecUtil (part of pySecDec) for unit testing.

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


## Examples

Several examples of pySecDec usage can be found in the source
code repository (`examples` directory). Here are some of them:

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

Note that the examples are not installed via the normal installation
procedure. They can only be found in the source code repository.

Development
===========

During development instead of full installation it is more
convenient to work directly from the checked out repository. To
make this work, one needs to install the Python dependencies,
build the contributed software, and setup PYTHONPATH to point
to the sources; do this by running

    $ make dependencies
    $ make build
    $ export PYTHONPATH=$PWD

Note that `make dependencies` will install the dependency
libraries using `python3` binary by default. A different Python
binary can be chosen for this by setting the `PYTHON` argument
to it, for example:

    $ make dependencies PYTHON=python3.10
    $ make build PYTHON=python3.10

The ``Makefile`` in the package's root directory also implements
other common development tasks.  You can list all available
targets with the command

    $ make help

`pySecDec` comes with a self test suite written in the `python
unittest` framework. The most convenient way to run all test is
using [pytest]. If pytest is installed (as it would be after `make
dependencies`), type

    $ make check

in the source repository to run all tests. Developers should write test cases for
ALL functions they implement and make sure that ALL tests pass before uploading a
commit.

[pytest]: https://docs.pytest.org/

## Building the Documentation

To build the documentation of `pySecDec`, you need [Sphinx]. If
Sphinx is installed (as it would be after `make dependencies`),
the command

    $ make doc

generates the documentation in `html` and in `pdf` format. Developers
should inline-document python functions and also keep the C++
part up to date. To generate HTML and PDF separately, use

    $ make doc-html

and

    $ make doc-pdf

[sphinx]: http://www.sphinx-doc.org

Building the documentation in pdf format requires an up-to-date
installation of a latex implementation. If you get an error about
missing `.sty` file, do the following:

 1. If you are an administrator on your computer, try to install
    the missing latex packages with your favorite package manager.
    The [TeXLive] or [MiKTeX] implementations should contain all
    required packages.

 2. If you are not an administrator, first get the missing
    packages, e.g. from [CTAN]. Collect the missing files in
    one or more directories and set the environment variable
    TEXINPUTS to these directories in the same way as the PATH
    variable is typically set.

[ctan]: http://www.ctan.org/
[miktex]: https://miktex.org/
[texlive]: http://tug.org/texlive/

## Making a PyPI release

To create a release on [PyPI], first increment the version number
of pySecDec in `pyproject.toml`, and update `CHANGELOG.md`. Then
make sure you have a recent [Singularity] installed, and run

    ./build-release.sh

This will create a source distribution file (`dist/*.tar.gz`)
and several prebuild distribution files (`dist/*.whl`); it will
also print the instructions on how to double-check the files
with [auditwheel] and upload them to PyPI using [Twine].

Note that an account with PyPI is needed for this. Please see
the [packaging python projects] guide for details.

The reason for using Singularity here is because prebuilt
distributions must be built inside one of the [manylinux] Docker
images: this guarantees that the prebuilt libraries and programs
will work on a wide range of Linux systems. If it wasn't for
this requirement, then making a distribution would be as simple
as running

    $ make dist

and then uploading the distfiles in `dist/`:

    $ twine upload dist/pySecDec-<version>.tar.gz
    $ twine upload dist/pySecDec-<version>-<tag>.whl

[packaging python projects]: https://packaging.python.org/tutorials/packaging-projects/
[twine]: https://pypi.org/project/twine/
[manylinux]: https://github.com/pypa/manylinux
[auditwheel]: https://pypi.org/project/auditwheel/
[singularity]: https://github.com/sylabs/singularity
