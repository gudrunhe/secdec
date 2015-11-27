"""The sector decomposition routines"""

from __future__ import print_function
from .polynomial import Polynomial, PolynomialProduct
from .sector import Sector
from .misc import powerset
import subprocess, shutil, os, re, numpy as np

# TODO: Split this file into two files - one for the iterative method and one for the geomethod.

# ********************** iterative decomposition **********************

def primary_decomposition_polynomial(polynomial):
    r'''
    Perform the primary decomposition on a single polynomial.

    .. seealso::
        :func:`.primary_decomposition`

    :param polynomial:
        :class:`.Polynomial`;
        The polynomial to eliminate the Dirac delta from.
    '''
    primary_sectors = []

    # get number of Feynman parameters
    N = polynomial.expolist.shape[1]

    coeffs = polynomial.coeffs
    expolist = polynomial.expolist

    for i in range(N):
        # "pinch" (delete) Feynman parameter `i`
        #   => this is equivalent to setting the exponent to zero
        #     => that is however equivalent to setting the parameter to one
        expolist_i = np.delete(expolist,i,axis=1)
        primary_sectors.append(Polynomial(expolist_i, coeffs))

    return primary_sectors

def primary_decomposition(sector):
    r'''
    Perform the primary decomposition as described in
    chapter 3.2 (part I) of arXiv:0803.4177v2.
    Return a list of :class:`.sector.Sector` - the primary
    sectors.
    For `N` Feynman parameters, there are `N` primary
    sectors where the `i`-th Feynman parameter is set to
    `1` in sector `i`.

    .. seealso::
        :func:`.primary_decomposition_polynomial`


    :param sector:
        :class:`.sector.Sector`;
        The container holding the polynomials (typically
        :math:`U` and :math:`F`) to eliminate the Dirac
        delta from.

    '''
    # get number of Feynman parameters
    N = sector.number_of_variables

    # Must perform `primary_decomposition_polynomial` for each element
    # of `sector.other`, `[sector.Jacobian]`, and each of the two factors
    # of `sector.cast`.
    primary_sector_polynomials_other = [primary_decomposition_polynomial(poly) for poly in sector.other]
    primary_sector_polynomials_Jacobian = primary_decomposition_polynomial(sector.Jacobian)
    primary_sector_polynomials_cast_factor0 = [primary_decomposition_polynomial(polyprod.factors[0]) for polyprod in sector.cast]
    primary_sector_polynomials_cast_factor1 = [primary_decomposition_polynomial(polyprod.factors[1]) for polyprod in sector.cast]

    # Collect the primary decomposed polynomials back into the `Sector` container class
    primary_sectors = [
        Sector(
                [PolynomialProduct(cast0[sector_index],cast1[sector_index]) for cast0,cast1 in zip(primary_sector_polynomials_cast_factor0, primary_sector_polynomials_cast_factor1)],
                [other[sector_index] for other in primary_sector_polynomials_other],
                primary_sector_polynomials_Jacobian[sector_index]
        ) for sector_index in range(N)
    ]
    return primary_sectors

class EndOfDecomposition(Exception):
    '''
    This exception is raised if the function
    :func:`.iteration_step` is called although
    the sector is already in standard form.

    '''
    pass

def refactorize(polyprod, parameter=None):
    '''
    In a :class:`.polynomial.PolynomialProduct` of
    the form `<monomial> * <polynomial>`, check if
    a parameter in `<polynomial>` can be shifted to
    the `<monomial>`.
    If possible, modify `polyprod` accordingly.

    :param polyprod:
        :class:`.polynomial.PolynomialProduct` of the
        form <monomial> * <polynomial>`;
        The product to refactorize.

    :param parameter:
        integer, optional;
        Check only the parameter with this index.
        If not provided, all parameters are checked.

    '''
    if parameter is None:
        for i in range(polyprod.number_of_variables):
            refactorize(polyprod, i)
        return

    expolist_mono = polyprod.factors[0].expolist
    expolist_poly = polyprod.factors[1].expolist

    factorizable_power = expolist_poly[:,parameter].min()

    expolist_mono[:,parameter] += factorizable_power
    expolist_poly[:,parameter] -= factorizable_power

def remap_parameters(singular_parameters, Jacobian, *polynomials):
    r'''
    Remap the Feynman parameters according to eq. (16) of
    arXiv:0803.4177v2. The parameter whose index comes first
    in `singular_parameters` is kept fix.

    The remapping is done in place; i.e. the `polynomials` are
    **NOT** copied.

    :param singular_parameters:
        list of integers;
        The indices :math:`\alpha_r` such that at least one
        of `polynomials` becomes zero if all
        :math:`t_{\alpha_r} \rightarrow 0`.

    :param Jacobian:
        :class:`.Polynomial` with one term and no coefficients;
        The Jacobian determinant is multiplied to this polynomial.

    :param polynomials:
        abritrarily many instances of :class:`.Polynomial` where
        all of these have an equal number of variables;
        The polynomials of Feynman parameters to be remapped. These
        are typically :math:`F` and :math:`U`.

    
    Example:

    .. code-block:: python

        remap_parameters([1,2], Jacobian, F, U)

    '''
    assert len(Jacobian.coeffs) == 1, "`Jacobian` must be a monomial."
    assert polynomials, "No polynomial for modification passed"

    num_parameters = polynomials[0].expolist.shape[1]

    for poly in polynomials:
        assert num_parameters == poly.expolist.shape[1], 'All polynomials must have the same number of variables'

    for poly in polynomials:
        for param in singular_parameters[1:]:
            poly.expolist[:,singular_parameters[0]] += poly.expolist[:,param] # This modifies in place!

    Jacobian.expolist[:,singular_parameters[0]] += len(singular_parameters) - 1

def iteration_step(sector):
    '''
    Run a single step of the iterative sector decomposition as described
    in chapter 3.2 (part II) of arXiv:0803.4177v2.
    Return an iterator of :class:`.Sector` - the arising subsectors.

    :param sector:
        :class:`.sector.Sector`;
        The sector to be decomposed.

    '''
    def get_poly_to_transform(sector):
        '''
        Return a :class:`PolynomialProduct` in `sector.cast`
        that is not in the desired form
        `<monomial> * <const + ...>` yet.
        Raise `EndOfDecomposition` if the desired form is
        reached.

        '''
        for polyprod in sector.cast:
            if not polyprod.factors[1].has_constant_term():
                return polyprod
        # Control only reaches this point if the desired form is
        # already reached.
        raise EndOfDecomposition()

    N = sector.number_of_variables

    # find a suitable transformation for a polynomial to be cast
    singular_set_found = False
    polyprod = get_poly_to_transform(sector)
    mono = polyprod.factors[0]
    poly = polyprod.factors[1]
    for singular_set in powerset(range(N),exclude_empty=True):
        if poly.becomes_zero_for(singular_set):
            singular_set_found = True
            break
    assert singular_set_found

    # We have to generate a subsector for each Feynman parameter
    # that appears in `singular_set`.
    # In order to comply with `remap_parameters`, we create
    # `len(singular_set)` copies of `singular_set` such that
    # each Feynman parameter is in the first position exactly
    # once.
    subsector_defining_singular_sets = [list(singular_set) for item in singular_set]
    for i,item in enumerate(subsector_defining_singular_sets):
        # swap the leading and the i-th item
        item[0],item[i] = item[i],item[0]

    # Call `remap_parameters` for each arising subsector.
    for singular_set in subsector_defining_singular_sets:
        subsector = sector.copy()
        remap_parameters(singular_set, subsector.Jacobian, *([polyprod.factors[1] for polyprod in subsector.cast] + subsector.other))
        for polyprod in subsector.cast:
            refactorize(polyprod,singular_set[0])
        yield subsector

def iterative_decomposition(sector):
    '''
    Run the iterative sector decomposition as described
    in chapter 3.2 (part II) of arXiv:0803.4177v2.
    Return an iterator of :class:`.Sector` - the arising subsectors.

    :param sector:
        :class:`.sector.Sector`;
        The sector to be decomposed.

    '''
    try:
        new_subsectors = iteration_step(sector)
        for subsector in new_subsectors:
            for deeper_subsector in iterative_decomposition(subsector):
                yield deeper_subsector
    except EndOfDecomposition:
        yield sector

# ********************** geometric decomposition **********************

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

    def complete_representation(self, normaliz='normaliz', workdir='normaliz_tmp', keep_workdir=False, verbose=False):
        '''
        Transform the vertex representation of a polytope
        to the facet representation or the other way round.
        Remove surplus entries in ``self.facets`` or
        ``self.vertices``.

        .. note::
            This function calls the command line executable of
            `normaliz <http://www.home.uni-osnabrueck.de/wbruns/
            normaliz/>`_..

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

        :param verbose:
            bool;
            Whether or not to print the command line messages
            of `normaliz`.

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
            normaliz_args = ['--ext', '--cst'] # create the files 'normaliz.ext' (vertices) and 'normaliz.cst' (facets)
            command_line_command = [normaliz] + normaliz_args + [run_card_filename]

            # dump run card to file
            with open(os.path.join(workdir, run_card_filename),'w') as f:
                f.write(run_card_as_str)

            if verbose:
                print('Normaliz run card (file "%s"):\n' % run_card_filename)
                print('-----------------------------------')
                print(run_card_as_str)
                print('-----------------------------------')
                print()
                print('running "%s" ...\n' % ' '.join(command_line_command))

            # run normaliz
            #    stdout=None and stderr=None --> print all normaliz output to the screen
            #    stdout=subprocess.PIPE and stderr=subprocess.PIPE --> redirect and discard output
            #    subprocess.check_call --> run normaliz, block until it finishes and raise error on nonzero exit status
            subprocess.check_call(command_line_command, stdout=None if verbose else subprocess.PIPE, stderr=None if verbose else subprocess.PIPE, cwd=workdir)

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
