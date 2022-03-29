"""
Polytope
--------

The polytope class as required by :mod:`pySecDec.make_regions`
and :mod:`pySecDec.decomposition.geometric`.

"""

from .algebra import Polynomial
import os, shutil, subprocess, re, numpy as np

import pySecDecContrib

def convex_hull(*polynomials):
    '''
    Calculate the convex hull of the Minkowski
    sum of all polynomials in the input.
    The algorithm sets all coefficients to one first
    and then only keeps terms of the polynomial product
    that have coefficient 1.
    Return the list of these entries in the expolist
    of the product of all input polynomials.

    :param polynomials:
        abritrarily many instances of :class:`.Polynomial` where
        all of these have an equal number of variables;
        The polynomials to calculate the convex hull for.

    '''
    product = 1

    for poly in polynomials:
        # want to use multiplication operator --> need type `Polynomial`, not a subtype
        # set the coefficients to one first --> terms with coefficients unequal to one are definitely not part of the convex hull
        product *= Polynomial(poly.expolist, np.ones_like(poly.coeffs, dtype=int))

    return product.expolist[np.where( product.coeffs == 1 )]

def normaliz_runcard(data, keyword, dimension):
    '''
    Returns string containing normaliz input file.

    :param data:
        two dimensional array;
        Input data.

    :param keyword:
        string;
        Keyword specifying the type of input data.
        For options see normaliz documentation.

    :param dimension:
        integer;
        Dimension of the ambient space.

    '''
    old_np_printoptions = np.get_printoptions()
    try:
        np.set_printoptions(threshold=np.inf)
        run_card_as_str  = 'amb_space ' + str(dimension) + '\n'
        run_card_as_str += keyword + ' ' + str(data.shape[0]) + '\n'
        run_card_as_str += str(data).replace('[','').replace(']','').replace('\n ','\n')
        return run_card_as_str
    finally:
        np.set_printoptions(**old_np_printoptions)

def read_normaliz_file(filepath, nofoutputs = 1):
    '''
    Read normaliz output.

    :param filepath:
        string;
        Normaliz output file to be read.

    :param nofoutputs:
        integer;
        Number of different types of output
        in the file to be read in.

    '''
    with open(filepath, 'r') as f:
        output = []
        for i in range(nofoutputs):
            array_as_str = ''
            # the first two lines contain the array dimensions
            shape = int(f.readline()), int(f.readline())
            for j in range(shape[0]):
                array_as_str += f.readline()
            # skip output type line
            f.readline()

            if shape[0] == 0:
                output.append(np.array([]))
            else:
                output.append(np.fromstring(array_as_str, sep=' ', dtype=int).reshape(shape))

    return output

def run_normaliz(normaliz=None, workdir='normaliz_tmp', run_card_filename='normaliz.in', normaliz_args=[]):
    '''
    Run normaliz.

    :param normaliz:
        string;
        The shell command to run `normaliz`.
        Default: use `normaliz` from pySecDecContrib

    :param workdir:
        string;
        The directory for the communication with `normaliz`.
        A directory with the specified name will be created
        in the current working directory. If the specified
        directory name already exists, an :class:`OSError`
        is raised.

        .. note::
            The communication with `normaliz` is done via
            files.

    :param run_card_filename:
        string;
        File name of normaliz input file.

    :param normaliz_args:
        list of strings;
        Normaliz command line arguments.

    '''
    if normaliz is None:
        normaliz = os.path.join(pySecDecContrib.dirname, 'bin', 'normaliz')
    command_line_command = [normaliz] + normaliz_args + [run_card_filename]

    # write additional information
    with open(os.path.join(workdir, 'run_info'),'w') as infofile:
        infofile.write('Normaliz run card: "%s"\n' % run_card_filename)
        infofile.write('running "%s" ...\n' % ' '.join(command_line_command))

    # redirect normaliz stdout
    with open(os.path.join(workdir, 'stdout'),'w') as stdout:
        # redirect normaliz stderr
        with open(os.path.join(workdir, 'stderr'),'w') as stderr:
            # run normaliz
            #    subprocess.check_call --> run normaliz, block until it finishes and raise error on nonzero exit status
            try:
                subprocess.check_call(command_line_command, stdout=stdout, stderr=stderr, cwd=workdir)
            except OSError as error:
                if normaliz not in str(error):
                    error.filename = normaliz
                raise

def triangulate(cone, normaliz=None, workdir='normaliz_tmp', keep_workdir=False, switch_representation=False):
    '''
    Split a cone into simplicial cones; i.e.
    cones defined by exactly :math:`D` rays
    where :math:`D` is the dimensionality.

    .. note::
        This function calls the command line executable of
        `normaliz` [BIR]_. See :ref:`installation_normaliz`
        for installation and a list of tested versions.


    :param cone:
        two dimensional array;
        The defining rays of the cone.

    :param normaliz:
        string;
        The shell command to run `normaliz`.
        Default: use `normaliz` from pySecDecContrib

    :param workdir:
        string;
        The directory for the communication with `normaliz`.
        A directory with the specified name will be created
        in the current working directory. If the specified
        directory name already exists, an :class:`OSError`
        is raised.

        .. note::
            The communication with `normaliz` is done via
            files.

    :param keep_workdir:
        bool;
        Whether or not to delete the `workdir` after execution.

    :param switch_representation:
        bool;
        Whether or not to switch between facet and vertex/ray
        representation.

    '''
    cone = np.asarray(cone)
    # basic consistency checks
    assert len(cone.shape) == 2, '`cone` must be two dimensional'
    assert cone.shape[0] >= cone.shape[1], 'Must at least have as many rays as the dimensionality'

    if cone.shape[0] == cone.shape[1] and not switch_representation:
        raise ValueError("`cone` is simplicial already")

    os.mkdir(workdir)
    try:
        if switch_representation == False:
            run_card_as_str = normaliz_runcard(cone, 'cone', cone.shape[1])
        else:
            run_card_as_str = normaliz_runcard(cone, 'inequalities', cone.shape[1])

        run_card_file_prefix = 'normaliz'
        run_card_file_suffix = '.in'
        run_card_filename = run_card_file_prefix + run_card_file_suffix

        # dump run card to file
        with open(os.path.join(workdir, run_card_filename),'w') as f:
            f.write(run_card_as_str)

        run_normaliz(normaliz=normaliz, workdir=workdir, run_card_filename=run_card_filename, normaliz_args=['-T', '--verbose']) # create the triangulation

        original_cone = read_normaliz_file(os.path.join(workdir,run_card_file_prefix + '.tgn'), 1)[0]

        # the triangulation is given as indices of `original_cone`
        if np.array_equal(original_cone, np.array([])):
            simplicial_cones_indices = []
        else:
            # `[:,:-1]` to delete the last column (last column are the determiants)
            # `-1` normaliz starts counting at `1` while python starts at `0`
            simplicial_cones_indices = read_normaliz_file(os.path.join(workdir,run_card_file_prefix + '.tri'), 1)[0][:,:-1] - 1

        return original_cone[simplicial_cones_indices]

    finally:
        if not keep_workdir:
            shutil.rmtree(workdir)

class Polytope(object):
    r'''
    Representation of a polytope defined by either its
    vertices or its facets. Call
    :meth:`.complete_representation` to translate from
    one to the other representation.

    :param vertices:
        two dimensional array;
        The polytope in vertex representation. Each
        row is interpreted as one vertex.

    :param facets:
        two dimensional array;
        The polytope in facet representation. Each
        row represents one facet :math:`F`.
        A row in `facets` is interpreted as one normal
        vector :math:`n_F` with additionally the constant
        :math:`a_F` in the last column.
        The points :math:`v` of the polytope obey

        .. math::
            \bigcap_F \left( {\langle n_F, v \rangle} + a_F \right) \ge 0

    :param equations:
        two dimensional array;
        Equations defining the hyperplanes the
        polytope is contained in. Only non-empty for
        polytopes that are not full-dimensional.

    '''
    def __init__(self, vertices=None, facets=None):
        if vertices is None and facets is None or vertices is not None and facets is not None:
            raise TypeError("Pass either `vertices` or `facets` but not both.")

        # copy and convert to numpy array
        if vertices is None:
            self.vertices = None
        else:
            self.vertices = np.array(vertices)

        if facets is None:
            self.facets = None
        else:
            self.facets = np.array(facets)

        self.equations = None

    def complete_representation(self, normaliz=None, workdir='normaliz_tmp', keep_workdir=False):
        '''
        Transform the vertex representation of a polytope
        to the facet representation or the other way round.
        Remove surplus entries in ``self.facets`` or
        ``self.vertices``.

        .. note::
            This function calls the command line executable of
            `normaliz` [BIR]_.

        :param normaliz:
            string;
            The shell command to run `normaliz`.
            Default: use `normaliz` from pySecDecContrib

        :param workdir:
            string;
            The directory for the communication with `normaliz`.
            A directory with the specified name will be created
            in the current working directory. If the specified
            directory name already exists, an :class:`OSError`
            is raised.

            .. note::
                The communication with `normaliz` is done via
                files.

        :param keep_workdir:
            bool;
            Whether or not to delete the `workdir` after execution.

        '''
        os.mkdir(workdir)
        try:
            if self.facets is None:
                run_card_as_str = normaliz_runcard(self.vertices, 'polytope', self.vertices.shape[1]+1)
            elif self.vertices is None:
                run_card_as_str = normaliz_runcard(self.facets, 'inequalities', self.facets.shape[1])
            else:
                raise ValueError('Both representations (facet and vertex) are already calculated')

            run_card_file_prefix = 'normaliz'
            run_card_file_suffix = '.in'
            run_card_filename = run_card_file_prefix + run_card_file_suffix

            # dump run card to file
            with open(os.path.join(workdir, run_card_filename),'w') as f:
                f.write(run_card_as_str)

            run_normaliz(normaliz=normaliz, workdir=workdir, run_card_filename=run_card_filename, normaliz_args=['--ext', '--cst', '-s', '--verbose']) # create the files 'normaliz.ext' (vertices) and 'normaliz.cst' (facets)

            # read file output from normaliz
            self.vertices = read_normaliz_file(os.path.join(workdir,run_card_file_prefix + '.ext'), 1)[0]
            self.facets, self.equations = read_normaliz_file(os.path.join(workdir,run_card_file_prefix + '.cst'), 2)

            # discard the last column that only consists of ones
            self.vertices = self.vertices[:,:-1]

        finally:
            if not keep_workdir:
                shutil.rmtree(workdir)

    def vertex_incidence_lists(self):
        '''
        Return for each vertex the list of facets it
        lies in (as dictonary). The keys of the output
        dictonary are the vertices while the values
        are the indices of the facets in ``self.facets``.

        '''
        assert self.vertices is not None and self.facets is not None, 'Run `complete_representation` first'
        # insert vertices in \bigcap_F \left( {\langle n_F, v \rangle} + a_F \right) \ge 0
        incidence_array = np.einsum('ij,kj', self.vertices, self.facets[:,:-1]) + self.facets[:,-1] == 0
        outdict = {}
        for i, vertex in enumerate(self.vertices):
            outdict[tuple(vertex)] = np.where(incidence_array[i])[0]
        return outdict
