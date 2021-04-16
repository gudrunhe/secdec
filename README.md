[![Documentation Status](https://readthedocs.org/projects/secdec/badge/?version=latest)](http://secdec.readthedocs.io/en/latest/?badge=latest)

Downloading the Program and Installation
========================================

pySecDec should run fine with both, python 2.7 and python 3
on unix-like systems.

Before you install pySecDec, make sure that you have
recent versions of `numpy` (http://www.numpy.org/) and
`sympy` (http://www.sympy.org/) installed. The version of `sympy`
should be 0.7.6 or higher, the version of `numpy` should be 1.6 or higher.

To install pySecDec download and upack the tarball from https://github.com/gudrunhe/secdec/releases.
The tarball contains a distribution of pySecDec and the additional dependencies
listed below. Typing

    $ make

should build all redistributed packages and display two commands
to be added to your .bashrc or .profile.

The geometric method and Normaliz
---------------------------------

If you want to use the geometric decomposition methods
(decomposition_method='geometric or 'geometric_ku'),
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

Additional dependencies for generated c++ packages
--------------------------------------------------

The intended main usage of pySecDec is to make it write c++ packages using the functions
`pySecDec.code_writer.make_package` and `pySecDec.loop_integral.loop_package`.
In order to build these c++ packages, the following additional non-python-based libraries
and programs are required:

 * CUBA (http://www.feynarts.de/cuba/)
 * FORM (http://www.nikhef.nl/~form/)
 * SecDecUtil (part of pySecDec), depends on:

   * catch (https://github.com/philsquared/Catch)
   * gsl (http://www.gnu.org/software/gsl/)

The functions `pySecDec.code_writer.make_package` and `pySecDec.loop_integral.loop_package`
can use the external program `nauty`
 to find all sector symmetries and therefore reduce the number of sectors:

 * NAUTY (http://pallini.di.uniroma1.it/)
[B. D. McKay and A. Piperno, Practical graph isomorphism, II, 2014, Journal of Symbolic Computation, 60, 94-112,
doi:10.1016/j.jsc.2013.09.003]

These packages are redistributed with the pySecDec tarball; i.e. you don't have to install
any of them yourself.


Basic Usage
===========

After installation, you should have a folder `examples` in your main pySecDec directory.
The full list of examples is given below.

A simple example of the evaluation of a loop integral with pySecDec is `box1L`.
This example computes a one-loop box with one off-shell leg (with
off-shellness s1) and one internal massive line (with mass squared msq).

To run the example change to the `box1L` directory and run the commands:

    $ python generate_box1L.py
    $ make -C box1L
    $ python integrate_box1L.py

The file generate_box1L.py defines the loop integral and calls pySecDec to perform the sector decomposition. When run it produces the directory `box1L` which contains the code required to numerically evaluate the integral. The make command builds this code and produces a library. The file integrate_box1L.py loads the integral library and evaluates the integral for a specified numerical point.

This will print the result of the integral evaluated with Mandelstam invariants s=4.0, t=-0.75 and s1=1.25, msq=1.0:

    eps^-2: -0.142868356275422825 - 1.63596224151119965e-6*I +/- ( 0.00118022544307414272 + 0.000210769456586696187*I )
    eps^-1: 0.639405625715768089 + 1.34277036689902802e-6*I +/- ( 0.00650722394065588166 + 0.000971496627153705891*I )
    eps^0 : -0.425514350373418893 + 1.86892487760861536*I +/- ( 0.00706834403694714484 + 0.0186497890361357298*I )


List of examples
================

 * easy:             a simple parametric integral, described in the manual in Section 2.1.
 * box1L:            a simple 1-loop, 4-point, 4-propagator integral, described in the manual Section 2.2.
 * triangle2L:   a 2-loop, 3-point, 6-propagator diagram, also containing massive propagators.
 * box2L_numerator: a massless planar on-shell 2-loop, 4-point, 7-propagator box with a numerator, either defined as an inverse propagator box2L_invprop.py or in terms of contracted Lorentz vectors box2L_contracted_tensor.py.
 * triangle3L:   a 2-loop, 3-point, 7-propagator integral, demonstrates that the symmetry finder can significantly reduce the number of sectors.
 * elliptic2L_euclidean:    an integral known to contain elliptic functions, evaluated at a Euclidean phase-space point.
 * elliptic2L_physical:     an integral known to contain elliptic functions, evaluated at a physical phase-space point.
 * triangle2L_split:    a 2-loop, 3-point, 6-propagator integral without a Euclidean region due to special kinematics.
 * hypergeo5F4:             a general dimensionally regulated parameter integral, corresponding to a Hypergeometric function 5F4.
 * 4photon1L_amplitude:     calculation of the 4-photon amplitude, showing how to use pySecDec as an integral library in a larger context.
 * two_regulators:  an integral involving poles in two different regulators.
 * userdefined_cpp:        a collection of examples demonstrating how to combine polynomials to be decomposed with other user-defined functions

Development
===========

The ``Makefile`` in the package's
root directory implements common development tasks.
You can list all available targets with the command

    $ make help

`pySecDec` comes with a self test suite written in the `python unittest` framework.
The most convenient way to run all test is using `nose` (http://nose.readthedocs.org).
If `nose` is installed, type

    $ make check

in the source repository to run all tests. Developers should write test cases for
ALL functions they implement and make sure that ALL tests pass before uploading a
commit.

Building the Documentation
--------------------------

To build the documentation of `pySecDec`, you need `sphinx` (http://www.sphinx-doc.org).
If `sphinx` is installed, the command

    $ make doc

generates the documentaion in `html` and in `pdf` format. Developers should inline-document
python functions and also keep the c++ part up to date.

Building the documentaion in pdf format requires an up-to-date installation of a latex
implementation. If you get an error about missing ".sty" file, do the following:

 1. If you are an administrator on your computer, try to install the missing latex packages
    with your favorite package manager. The MiKTeX or TeXLive implementations should contain
    all required packages.

 2. If you are not an administrator, first get the missing packages, e.g. from
    "http://www.ctan.org/". Collect the missing files in one or more directories
    and set the environment variable TEXINPUTS to these directories in the same
    way as the PATH variable is typically set.

Making a tarball
----------------

A distributable tarball can be created by the command

    $ make dist

When publishing a new release, make sure that you increase the Version number of `pySecDec`
and/or `SecDecUtil`. You should also describe the main changes compared to the previous release
in the `ChangeLog`. Further, make sure that you do not have any generated code in your
`examples` directory: It will go into the tarball if present!

"make dist" first runs "make clean" to make sure that no temporary or outdated files go into
the distribution. Then it runs builds the source distribution of the python package as with
"python setup.py sdist", the documentation as with "make doc", and it runs "make dist" in the
"util" package. Note that the "util" is an autotools package and you must have the GNU autotools
installed on your system.

Installation for developers
---------------------------

This section describes how to install pySecDec and easily update to the commit that is currently
checked out in the git repository. The following steps should be seen as a guideline to achieve this:

1) Clone the git repository.

2) Create the tarball by running "make dist". If that does not work out of the box, make sure that you
   have recent versions of the following packages installed:
   * python (https://www.python.org/)
   * Sphinx (http://www.sphinx-doc.org/)
   * The full TeX Live distribution (http://tug.org/texlive/)
     Note that there are slimmed distributions of TeXLive that typically result in errors about missing files.
     Follow the instructions in the section "Building the Documentation" in that case.
   * The GNU autotools:
     * Autoconf(https://www.gnu.org/software/autoconf/autoconf.html)
     * Automake(https://www.gnu.org/software/automake/)
     * Libtool(https://www.gnu.org/software/libtool/)

3) Unpack the tarball to a location OUTSIDE of the repository.

4) Run "make" in the directory with the unpacked tarball.

5) The success message tells you two environment variables to be set in the ".profile"/".bashrc". Set SECDEC_CONTRIB
   as shown in the message, but set PYTHONPATH to the root directory of the git repository. With that, python will
   always load the version of pySecDec that is currently checked out in the git directory.

6) Open a new shell and make sure that the the environment variables PYTHONPATH and SECDEC_CONTRIB are set by typing
   "echo ${SECDEC_CONTRIB}" and "echo ${PYTHONPATH}".

7) In the shell with the new variables set, cd to the git repository and then into "util".
   Run "./configure --prefix=$SECDEC_CONTRIB && make install" there.

After following the steps above, you must run "make install" in the "util" directory whenever the secdecutil is updated.

