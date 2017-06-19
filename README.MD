Development
-----------

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

Documentation
-------------

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

1) Clone the git repository from hepforge (ask another developer for details about that)

2) Create the tarball by running "make dist". If that does not work out of the box, make sure that you
   have recent versions of the following packages installed:
   * python (https://www.python.org/)
   * Sphinx (http://www.sphinx-doc.org/)
   * The full TeX Live distribution (http://tug.org/texlive/)
     Note that there are slimmed distributions of TeXLive that typically result in errors about missing files.
     Follow the instructions in the section "Documentation" in that case.
   * The GNU autotools:
     * Autoconf(https://www.gnu.org/software/autoconf/autoconf.html)
     * Automake(https://www.gnu.org/software/automake/)
     * Libtool(https://www.gnu.org/software/libtool/)

3) Unpack the tarball to a location OUTSIDE of the repository.

4) Run "make" in the directory with the unpacked tarball. Note that makes's "-j" option sometimes causes problems.
   If you don't see a success message after running "make -j<jmake>", rerun without the "-j" option.

5) The success message tells you two environment variables to be set in the ".profile"/".bashrc". Set SECDEC_CONTRIB
   as shown in the message, but set PYTHONPATH to the root directory of the git repository. With that, python will
   always load the version of pySecDec that is currently checked out in the git directory.

6) Open a new shell and make sure that the the environment variables PYTHONPATH and SECDEC_CONTRIB are set by typing
   "echo ${SECDEC_CONTRIB}" and "echo ${PYTHONPATH}".

7) In the shell with the new variables set, cd to the git repository and then into "util".
   Run "./configure --prefix=$SECDEC_CONTRIB && make install" there.

After following the steps above, you must run "make install" in the "util" directory whenever the secdecutil is updated.
