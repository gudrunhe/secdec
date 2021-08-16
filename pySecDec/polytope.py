"""
Polytope
--------

The polytope class as required by :mod:`pySecDec.make_regions`
and :mod:`pySecDec.decomposition.geometric`.

"""

from .algebra import Polynomial
import os, shutil, subprocess, re, numpy as np

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

def triangulate(cone, normaliz='normaliz', workdir='normaliz_tmp', keep_workdir=False, switch_representation=False):
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

    old_np_printoptions = np.get_printoptions()
    np.set_printoptions(threshold=np.inf)

    os.mkdir(workdir)
    try:
        # generate the normaliz run card
        run_card_as_str  = str(cone.shape[0]) + ' ' + str(cone.shape[1]) + '\n'
        run_card_as_str += str(cone).replace('[','').replace(']','').replace('\n ','\n')
        if switch_representation == False:
            run_card_as_str += '\nintegral_closure\n'
        else:
            run_card_as_str += '\ninequalities\n'

        run_card_file_prefix = 'normaliz'
        run_card_file_suffix = '.in'
        run_card_filename = run_card_file_prefix + run_card_file_suffix
        normaliz_args = ['-T', '--verbose'] # create the triangulation
        command_line_command = [normaliz] + normaliz_args + [run_card_filename]

        # dump run card to file
        with open(os.path.join(workdir, run_card_filename),'w') as f:
            f.write(run_card_as_str)

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

        # read normaliz output
        # normaliz reorders the rays and defines its ordering in "normaliz.tgn"
        with open(os.path.join(workdir, 'normaliz.tgn'),'r') as f:
            # the first two lines may contain the array dimensions
            try:
                shape = int(f.readline()), int(f.readline())
            except ValueError:
                # more than a single number in first two lines
                # --> file does not contain dimensions
                f.seek(0)
            original_cone = np.loadtxt(f, dtype=int, ndmin = 2)

        # the triangulation is given as indices of `original_cone`
        with open(os.path.join(workdir, 'normaliz.tri'),'r') as f:
            # the first two lines may contain the array dimensions
            try:
                shape = int(f.readline()), int(f.readline())
            except ValueError:
                # more than a single number in first two lines
                # --> file does not contain dimensions
                f.seek(0)

            # it is terminated by a line containing letters
            array_lines = []
            current_str = f.readline()
            while re.match(r'^[\-0-9 ]+$', current_str) is not None:
                array_lines.append(current_str)
                current_str = f.readline()
            array_as_str = ''.join(array_lines)

        # `[:,:-1]` to delete the last column (last column are the determiants)
        # `-1` normaliz starts counting at `1` while python starts at `0`
        simplicial_cones_indices = np.fromstring(array_as_str, sep=' ', dtype=int).reshape(len(array_lines),-1)[:,:-1] - 1

        return original_cone[simplicial_cones_indices]

    finally:
        np.set_printoptions(**old_np_printoptions)
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

    def complete_representation(self, normaliz='normaliz', workdir='normaliz_tmp', keep_workdir=False):
        '''
        Transform the vertex representation of a polytope
        to the facet representation or the other way round.
        Remove surplus entries in ``self.facets`` or
        ``self.vertices``.

        .. note::
            This function calls the command line executable of
            `normaliz` [BIR]_. See :ref:`installation_normaliz`
            for installation and a list of tested versions.

        :param normaliz:
            string;
            The shell command to run `normaliz`.

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
                run_card_as_str = self._make_run_card_vertices2facets()
            elif self.vertices is None:
                run_card_as_str = self._make_run_card_facets2vertices()
            else:
                raise ValueError('Both representations (facet and vertex) are already calculated')

            run_card_file_prefix = 'normaliz'
            run_card_file_suffix = '.in'
            run_card_filename = run_card_file_prefix + run_card_file_suffix
            normaliz_args = ['--ext', '--cst', '--verbose'] # create the files 'normaliz.ext' (vertices) and 'normaliz.cst' (facets)
            command_line_command = [normaliz] + normaliz_args + [run_card_filename]

            # dump run card to file
            with open(os.path.join(workdir, run_card_filename),'w') as f:
                f.write(run_card_as_str)

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

            # read file output from normaliz
            if self.facets is None: # reduced vertices are in cst file
                vertices = self._read_ext_file(os.path.join(workdir,run_card_file_prefix + '.ext'))
                self.facets = self._read_cst_file(os.path.join(workdir,run_card_file_prefix + '.cst'))
            elif self.vertices is None: # reduced facets are in cst file
                self.facets = self._read_ext_file(os.path.join(workdir,run_card_file_prefix + '.ext'))
                vertices = self._read_cst_file(os.path.join(workdir,run_card_file_prefix + '.cst'))
            else:
                raise RuntimeError()

            # discard the last column that only consists of ones
            self.vertices = vertices[:,:-1]

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

    def _make_run_card_vertices2facets(self):
        old_np_printoptions = np.get_printoptions()
        try:
            np.set_printoptions(threshold=np.inf)
            run_card_as_str  = str(self.vertices.shape[0]) + ' ' + str(self.vertices.shape[1]) + '\n\n'
            run_card_as_str += str(self.vertices).replace('[',' ').replace(']',' ')
            run_card_as_str += '\n\npolytope\n'
            return run_card_as_str
        finally:
            np.set_printoptions(**old_np_printoptions)

    def _make_run_card_facets2vertices(self):
        old_np_printoptions = np.get_printoptions()
        try:
            np.set_printoptions(threshold=np.inf)
            run_card_as_str  = str(self.facets.shape[0]) + ' ' + str(self.facets.shape[1]) + '\n'
            run_card_as_str += str(self.facets).replace('[','').replace(']','').replace('\n ','\n')
            run_card_as_str += 'integral_closure\n'
            return run_card_as_str
        finally:
            np.set_printoptions(**old_np_printoptions)

    def _read_cst_file(self, filepath):
        with open(filepath, 'r') as f:
            # the first two lines may contain the array dimensions
            try:
                shape = int(f.readline()), int(f.readline())
            except ValueError:
                # more than a single number in first two lines
                # --> file does not contain dimensions
                f.seek(0)
            # the reduced input comes next (space and newline separated) and is terminated by a line containing letters
            array_lines = []
            current_str = f.readline()
            while 'inequalities' not in current_str:
                array_lines.append(current_str)
                current_str = f.readline()
            array_as_str = ''.join(array_lines)
            # check for equations that constrain the polytope to hyperplane (happens for scaleless integrals)
            # if there are no equations the line is expected to be '0' or 'equations' depending on the normaliz version
            current_str= f.readline()

            if not ('equations' in current_str or re.sub(r'[\n\t\s]*', '', current_str) == '0'):
                raise NotImplementedError("Polytope is not full dimensional. Are you trying to compute a scaleless integral (which evaluates to zero)?")

            return np.fromstring(array_as_str, sep=' ', dtype=int).reshape(len(array_lines),-1)


    def _read_ext_file(self, filepath):
        with open(filepath, 'r') as f:
            # the first two lines may contain the array dimensions
            try:
                shape = int(f.readline()), int(f.readline())
            except ValueError:
                # more than a single number in first two lines
                # --> file does not contain dimensions
                f.seek(0)
            return np.loadtxt(f, dtype=int)
