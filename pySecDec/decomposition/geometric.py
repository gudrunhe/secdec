"""

The geometric sector decomposition routines.

"""

from .common import Sector, refactorize
from ..algebra import Polynomial, Product
import subprocess, shutil, os, re, itertools, numpy as np

# *********************** primary decomposition ***********************

def Cheng_Wu(sector, index=-1):
    '''
    Replace one Feynman parameter by one.
    This means integrating out the Dirac
    delta according to the Cheng-Wu theorem.

    :param sector:
        :class:`.Sector`;
        The container holding the polynomials (typically
        :math:`U` and :math:`F`) to eliminate the Dirac
        delta from.

    :param index:
        integer, optional;
        The index of the Feynman parameter to eliminate.
        Default: -1 (the last Feynman parameter)

    '''
    # do not remove the only parameter from the `polysymbols`
    remove = sector.number_of_variables != 1

    Jacobian = sector.Jacobian.replace(index, 1, remove)
    other = [poly.replace(index, 1, remove) for poly in sector.other]
    cast = [Product( *(product.factors[i].replace(index, 1, remove) for i in (0,1)) ) for product in sector.cast]
    return Sector(cast, other, Jacobian)

# ********************** geometric decomposition **********************
def generate_fan(*polynomials):
    '''
    Calculate the fan of the polynomials in the input. The rays of a
    cone are given by the exponent vectors after factoring out a monomial
    together with the standard basis vectors. Each choice of factored out
    monomials gives a different cone.
    Only full (:math:`N`-) dimensional cones in :math:`R^N_{\geq 0}` need to be
    considered.

    :param polynomials:
        abritrarily many instances of :class:`.Polynomial` where
        all of these have an equal number of variables;
        The polynomials to calculate the fan for.
    '''
    expolists = [poly.expolist for poly in polynomials]
    factors = itertools.product(*expolists)
    number_of_variables = polynomials[0].number_of_variables
    identity_polynomial = Polynomial(np.identity(number_of_variables, dtype=int), np.ones(number_of_variables, dtype=int), copy=False)
    fan = []

    for factor in factors:
        factor = np.array(factor)

        # reshape to use numpy's broadcasting
        cone = np.vstack( expolists[i] - factor[i].reshape(1,number_of_variables) for i in range(len(polynomials)) )

        # use `Polynomial` class to remove duplicates
        cone_poly = Polynomial(cone, np.ones(len(cone),dtype=int), copy=False)
        cone_poly += identity_polynomial # implicit simplify
        for i,hyperplane in enumerate(cone_poly.expolist):
            if (hyperplane > 0).all() or (hyperplane == 0).all():
                cone_poly.coeffs[i] = 0
        cone = cone_poly.simplify().expolist

        if (
                len(cone) >= number_of_variables and
                not (cone < 0).all(axis=1).any() and # if one hyperplane has only negative entries do not append to fan
                not any( (hyperplane==cone).all(axis=-1).any() for hyperplane in -cone ) # do not append cones that have `hyperplane` and `-hyperplane`
           ):
                fan.append(cone)

    return fan

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
            # the first two lines contain the array dimensions
            shape = int(f.readline()), int(f.readline())
            original_cone = np.fromfile(f, sep=' ', dtype=int).reshape(shape)

        # the triangulation is given as indices of `original_cone`
        with open(os.path.join(workdir, 'normaliz.tri'),'r') as f:
            # the first two lines contain the array dimensions
            shape = int(f.readline()), int(f.readline())

            # it is terminated by a line containing letters
            array_as_str = ''
            current_str = f.readline()
            while re.match(r'^[\-0-9 ]+$', current_str) is not None:
                array_as_str += current_str
                current_str = f.readline()

        # `[:,:-1]` to delete the last column (last column are the determiants)
        # `-1` normaliz starts counting at `1` while python starts at `0`
        simplicial_cones_indices = np.fromstring(array_as_str, sep=' ', dtype=int).reshape(shape)[:,:-1] - 1

        return original_cone[simplicial_cones_indices]

    finally:
        np.set_printoptions(**old_np_printoptions)
        if not keep_workdir:
            shutil.rmtree(workdir)

def transform_variables(polynomial, transformation, polysymbols='y'):
    r'''
    Transform the parameters :math:`x_i` of a
    :class:`pySecDec.algebra.Polynomial`,

    .. math::
        x_i \rightarrow \prod_j x_j^{T_{ij}}

    , where :math:`T_{ij}` is the transformation matrix.

    :param polynomial:
        :class:`pySecDec.algebra.Polynomial`;
        The polynomial to transform the variables in.

    :param transformation:
        two dimensional array;
        The transformation matrix :math:`T_{ij}`.

    :param polysymbols:
        string or iterable of strings;
        The symbols for the new variables. This argument
        is passed to the default constructor of
        :class:`pySecDec.algebra.Polynomial`.
        Refer to the documentation of
        :class:`pySecDec.algebra.Polynomial`
        for further details.

    '''
    new_expolist = polynomial.expolist.dot(transformation)
    number_of_new_variables = transformation.shape[-1]
    if isinstance(polysymbols, str):
        polysymbols = [polysymbols + str(i) for i in range(number_of_new_variables)]

    # keep the type (`polynomial` can have a subtype of `Polynomial`)
    outpoly = polynomial.copy()
    outpoly.expolist = new_expolist
    outpoly.polysymbols = polysymbols
    outpoly.number_of_variables = number_of_new_variables
    return outpoly

def geometric_decomposition(sector, indices=None, normaliz='normaliz', workdir='normaliz_tmp'):
    '''
    Run the sector decomposition using the geomethod
    as described in [BHJ+15]_.

    .. note::
        This function calls the command line executable of
        `normaliz` [BIR]_. See :ref:`installation_normaliz`
        for installation and a list of tested versions.

    :param sector:
        :class:`.Sector`;
        The sector to be decomposed.

    :param indices:
        list of integers or None;
        The indices of the parameters to be considered as
        integration variables. By default (``indices=None``),
        all parameters are considered as integration
        variables.

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

    '''
    original_sector = sector
    sector = original_sector.copy()

    if indices is None:
        indices = range(original_sector.number_of_variables)
    else:
        # remove parameters that are not in `indices`
        indices = list(indices)
        sector.number_of_variables = len(indices)
        sector.Jacobian.number_of_variables = len(indices)
        sector.Jacobian.expolist = sector.Jacobian.expolist[:,indices]
        sector.Jacobian.polysymbols = [sector.Jacobian.polysymbols[i] for i in indices]
        for product in sector.cast:
            for factor in product.factors:
                factor.number_of_variables = len(indices)
                factor.expolist = factor.expolist[:,indices]
                factor.polysymbols = [factor.polysymbols[i] for i in indices]
        for poly in sector.other:
            poly.number_of_variables = len(indices)
            poly.expolist = poly.expolist[:,indices]
            poly.polysymbols = [poly.polysymbols[i] for i in indices]

    dim = sector.number_of_variables

    polytope_vertices = convex_hull( *(product.factors[1] for product in sector.cast) )
    polytope = Polytope(vertices=polytope_vertices)
    polytope.complete_representation(normaliz, workdir)

    transformation = polytope.facets.T[:-1] # do not need offset term "a_F"
    incidence_lists = polytope.vertex_incidence_lists()

    # transform the variables for every polynomial of the sector
    sector.Jacobian = transform_variables(sector.Jacobian, transformation, sector.Jacobian.polysymbols)
    for i,product in enumerate(sector.cast):
        transformed_monomial = transform_variables(product.factors[0], transformation, product.factors[0].polysymbols)
        transformed_polynomial = transform_variables(product.factors[1], transformation, product.factors[1].polysymbols)
        sector.cast[i] = Product(transformed_monomial, transformed_polynomial)
    for i,polynomial in enumerate(sector.other):
        sector.other[i] = transform_variables(polynomial, transformation, polynomial.polysymbols)
    # this transformation produces an extra Jacobian factor
    # can multiply part encoded in the `expolist` here but the coefficient is specific for each subsector
    sector.Jacobian *= Polynomial([transformation.sum(axis=0) - 1], [1])

    def make_sector(cone_indices, cone):
        subsector = original_sector.copy()
        Jacobian_coeff = abs(np.linalg.det(cone))
        Jacobian_coeff_as_int = int(Jacobian_coeff + 0.5) # `Jacobian_coeff` is integral but numpy calculates it as float
        assert abs(Jacobian_coeff_as_int - Jacobian_coeff) < 1.0e-5 * abs(Jacobian_coeff)
        subsector.Jacobian *= Jacobian_coeff_as_int

        # set variables to one that are not in `cone_indices`
        number_of_variables = len(cone_indices) + original_sector.number_of_variables - len(indices)
        assert number_of_variables == original_sector.number_of_variables
        subsector.Jacobian.expolist[:,indices] = sector.Jacobian.expolist[:,cone_indices]
        for resulting_product, output_product in zip(sector.cast, subsector.cast):
            for j in range(2):
                output_product.factors[j].expolist[:,indices] = resulting_product.factors[j].expolist[:,cone_indices]
            refactorize(output_product)
        for resulting_polynomial, output_polynomial in zip(sector.other, subsector.other):
            output_polynomial.expolist[:,indices] = resulting_polynomial.expolist[:,cone_indices]

        return subsector


    for cone_indices in incidence_lists.values():
        cone = transformation[:,cone_indices].T

        # triangluate where neccessary
        if len(cone_indices) != dim:
            # assert len(cone) > dim # --> this check is done by `triangulate`
            triangular_cones = triangulate(cone, normaliz, workdir)

            assert len(triangular_cones.shape) == 3
            for i, triangular_cone in enumerate(triangular_cones):
                triangular_cone_indices = []
                for vector in triangular_cone:
                    # find the indices of the vectors defining the triangular cone
                    triangular_cone_indices.append(int(  np.where( (vector == transformation.T).all(axis=1) )[0]  ))
                yield make_sector(triangular_cone_indices, triangular_cone)

        else:
            yield make_sector(cone_indices, cone)

def geometric_decomposition_ku(sector, indices=None, normaliz='normaliz', workdir='normaliz_tmp'):
    '''
    Run the sector decomposition using the original geometric
    decomposition strategy by Kaneko and Ueda as described
    in [KU10]_.

    .. note::
        This function calls the command line executable of
        `normaliz` [BIR]_. See :ref:`installation_normaliz`
        for installation and a list of tested versions.

    :param sector:
        :class:`.Sector`;
        The sector to be decomposed.

    :param indices:
        list of integers or None;
        The indices of the parameters to be considered as
        integration variables. By default (``indices=None``),
        all parameters are considered as integration
        variables.

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

    '''
    original_sector = sector
    sector = original_sector.copy()

    if indices is None:
        indices = range(original_sector.number_of_variables)
    else:
        # remove parameters that are not in `indices`
        indices = list(indices)
        sector.number_of_variables = len(indices)
        sector.Jacobian.number_of_variables = len(indices)
        sector.Jacobian.expolist = sector.Jacobian.expolist[:,indices]
        sector.Jacobian.polysymbols = [sector.Jacobian.polysymbols[i] for i in indices]
        for product in sector.cast:
            for factor in product.factors:
                factor.number_of_variables = len(indices)
                factor.expolist = factor.expolist[:,indices]
                factor.polysymbols = [factor.polysymbols[i] for i in indices]
        for poly in sector.other:
            poly.number_of_variables = len(indices)
            poly.expolist = poly.expolist[:,indices]
            poly.polysymbols = [poly.polysymbols[i] for i in indices]

    def make_sector_ku(cone):
        subsector = original_sector.copy()
        transformation = np.identity(original_sector.number_of_variables, dtype = int)
        index_array = np.array(indices)
        transformation[index_array.reshape(-1,1),index_array] = cone

        Jacobian_coeff = abs(np.linalg.det(cone))
        Jacobian_coeff_as_int = int(Jacobian_coeff + 0.5) # `Jacobian_coeff` is integral but numpy calculates it as float
        assert abs(Jacobian_coeff_as_int - Jacobian_coeff) < 1.0e-5 * abs(Jacobian_coeff)

        subsector.Jacobian = Jacobian_coeff_as_int*transform_variables(subsector.Jacobian, transformation, subsector.Jacobian.polysymbols)

        for i,product in enumerate(subsector.cast):
            transformed_monomial = transform_variables(product.factors[0], transformation, product.factors[0].polysymbols)
            transformed_polynomial = transform_variables(product.factors[1], transformation, product.factors[1].polysymbols)
            subsector.cast[i] = Product(transformed_monomial, transformed_polynomial)
            refactorize(subsector.cast[i])
        for i,polynomial in enumerate(subsector.other):
            subsector.other[i] = transform_variables(polynomial, transformation, polynomial.polysymbols)
        # this transformation produces an extra Jacobian factor
        # can multiply part encoded in the `expolist` here but the coefficient is specific for each subsector
        subsector.Jacobian *= Polynomial([transformation.sum(axis=0) - 1], [1])

        return subsector

    fan = generate_fan( *(product.factors[1] for product in sector.cast) )
    for cone in fan:
        for dualcone in triangulate(cone, normaliz, workdir, switch_representation=True):
            # exclude lower dimensional cones
            if dualcone.shape[0] == cone.shape[1]:
                yield make_sector_ku(dualcone.T)
