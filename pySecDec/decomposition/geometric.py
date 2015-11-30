"""The geometric sector decomposition routines"""

import subprocess, shutil, os, re, numpy as np

# ********************** geometric decomposition **********************

def Cheng_Wu(sector, index=-1):
    '''
    Replace one Feynman parameter by one.
    This means integrating out the Dirac
    delta according to the Cheng-Wu theorem.

    :param sector:
        :class:`.sector.Sector`;
        The container holding the polynomials (typically
        :math:`U` and :math:`F`) to eliminate the Dirac
        delta from.

    :param index:
        integer, optional;
        The index of the Feynman parameter to eliminate.
        Default: -1 (the last Feynman parameter)

    '''
    raise NotImplementedError()
    # TODO implement and test

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
        poly = poly.copy()
        poly.coeffs = np.ones_like(poly.coeffs, dtype=int) # set coeffs to one
        product *= poly

    return product.expolist[np.where( product.coeffs == 1 )]

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
            `normaliz <http://www.home.uni-osnabrueck.de/wbruns/
            normaliz/>`_.
            It is designed for `normaliz` version 3.0.0

        :param normaliz:
            string or None;
            If this is None, return the `normaliz` run card
            as string.
            Otherwise, this specifies the shell command to
            run `normaliz`.

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
                infofile.write('Normaliz run card (file "%s"):\n' % run_card_filename)
                infofile.write('-----------------------------------\n')
                infofile.write(run_card_as_str)
                infofile.write('\n-----------------------------------\n')
                infofile.write('\n')
                infofile.write('running "%s" ...\n' % ' '.join(command_line_command))

            # redirect normaliz stdout
            with open(os.path.join(workdir, 'stdout'),'w') as stdout:
                # redirect normaliz stderr
                with open(os.path.join(workdir, 'stderr'),'w') as stderr:
                    # run normaliz
                    #    subprocess.check_call --> run normaliz, block until it finishes and raise error on nonzero exit status
                    subprocess.check_call(command_line_command, stdout=stdout, stderr=stderr, cwd=workdir)

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

    def _make_run_card_vertices2facets(self):
        run_card_as_str  = str(self.vertices.shape[0]) + ' ' + str(self.vertices.shape[1]) + '\n\n'
        run_card_as_str += str(self.vertices).replace('[',' ').replace(']',' ')
        run_card_as_str += '\n\npolytope\n'
        return run_card_as_str

    def _make_run_card_facets2vertices(self):
        run_card_as_str  = str(self.facets.shape[0]) + ' ' + str(self.facets.shape[1]) + '\n'
        run_card_as_str += str(self.facets).replace('[','').replace(']','').replace('\n ','\n')
        run_card_as_str += 'integral_closure\n'
        return run_card_as_str

    def _read_cst_file(self, filepath):
        with open(filepath, 'r') as f:
            # the first two lines contain the array dimensions
            shape = int(f.readline()), int(f.readline())

            # the reduced input comes next (space and newline separated) and is terminated by a line containing letters
            array_as_str = ''
            current_str = f.readline()
            while re.match(r'^[\-0-9 ]+$', current_str) is not None:
                array_as_str += current_str
                current_str = f.readline()

        return np.fromstring(array_as_str, sep=' ', dtype=int).reshape(shape)

    def _read_ext_file(self, filepath):
        with open(filepath, 'r') as f:
            # the first two lines contain the array dimensions
            shape = int(f.readline()), int(f.readline())
            return np.fromfile(f, sep=' ', dtype=int).reshape(shape)
