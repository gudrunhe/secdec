'''
This module collects routines that are used by
multiple decompition modules.

'''

from ..algebra import Polynomial, ExponentiatedPolynomial
from ..misc import argsort_ND_array, sympify_expression
import numpy as np
import sympy as sp
import subprocess, shutil, os

class Sector(object):
    '''
    Container class for sectors that arise during the
    sector decomposition.

    :param cast:
        iterable of :class:`.algebra.Product` or
        of :class:`.algebra.Polynomial`;
        The polynomials to be cast to the form
        `<monomial> * (const + ...)`

    :param other:
        iterable of :class:`.algebra.Polynomial`, optional;
        All variable transformations are applied to these
        polynomials but it is not attempted to achieve the
        form `<monomial> * (const + ...)`

    :param Jacobian:
        :class:`.algebra.Polynomial` with one term, optional;
        The Jacobian determinant of this sector. If not provided,
        the according unit monomial (1*x0^0*x1^0...) is assumed.

    '''
    def __init__(self, cast, other=[], Jacobian=None):
        cast  = list(cast)
        other = list(other)

        assert cast, "Must have at least one input polynomial"

        # The elements of `cast` may be of type `Polynomial` or `Product`
        try:
            assert len(cast[0].factors) == 2, "Every `Product` must have exactly two factors" # raises AttributeError if type is `Polynomial`
            poly = cast[0].factors[1]
        except AttributeError:
            poly = cast[0]
        self.number_of_variables = poly.expolist.shape[1]

        initial_monomial_factor = Polynomial([ [0]*self.number_of_variables ], [1], poly.polysymbols)

        if Jacobian is None:
            Jacobian = initial_monomial_factor

        for item in cast + other + [Jacobian]:
            # explanation see above
            try:
                assert len(item.factors) == 2, "Every `Product` must have exactly two factors" # raises AttributeError if type is `Polynomial`
                # if control reaches this point, assume that `item` has type `Product`
                poly = item.factors[1]
                assert len(item.factors[0].coeffs) == 1, 'The first factor of every `Product` must be a monomial'
            except AttributeError:
                poly = item
            # number of variables must match for all input polynomials
            assert self.number_of_variables == poly.expolist.shape[1], "The number of variables must be equal for all input polynomials"

        self.Jacobian = Jacobian.copy()

        self.other = [poly.copy() for poly in other]
        for poly in self.other:
            assert isinstance(poly, Polynomial), "All elements in `other` must be of type `Polynomial`"

        self.cast = []
        for item in cast:
            if hasattr(item, 'factors'): # expect type `Product`
                self.cast.append(item.copy())
            else: # expect type `Polynomial`
                self.cast.append(item.refactorize())

    def __repr__(self):
        return 'Sector(Jacobian=%s, cast=%s, other=%s)' % (self.Jacobian, self.cast, self.other)

    def __str__(self):
        return 'Sector:\nJacobian=%s\ncast=%s\nother=%s' % (self.Jacobian, self.cast, self.other)

    def copy(self):
        "Return a copy of a :class:`.Sector`."
        return Sector(self.cast, self.other, self.Jacobian)

def refactorize(polyprod, parameter=None):
    '''
    In a :class:`.algebra.Product` of
    the form `<monomial> * <polynomial>`, check if
    a parameter in `<polynomial>` can be shifted to
    the `<monomial>`.
    If possible, modify `polyprod` accordingly.

    :param polyprod:
        :class:`.algebra.Product` of the
        form <monomial> * <polynomial>`;
        The product to refactorize.

    :param parameter:
        integer, optional;
        Check only the parameter with this index.
        If not provided, all parameters are checked.

    '''
    expolist_mono = polyprod.factors[0].expolist
    expolist_poly = polyprod.factors[1].expolist

    if parameter is None:
        factorizable_powers = expolist_poly.min(axis=0)
        expolist_mono[:] += factorizable_powers
        expolist_poly[:] -= factorizable_powers
    else:
        factorizable_power = expolist_poly[:,parameter].min()
        expolist_mono[:,parameter] += factorizable_power
        expolist_poly[:,parameter] -= factorizable_power

# -------------------- finding symmetries ---------------------

def _sector2array(sector):
    '''
    Combine the `expolist`s and the `coeff`s
    of all :class:`.Polynomial`s in a
    :class:`.Sector` to two large arrays.
    Return the combined expolists and the
    combined coeffs.

    :param sector:
        :class:`.Sector`; The container of the
        :class:`.Polynomial`s to be combined.

    '''
    # process `Jacobian`
    combined_expolists = [sector.Jacobian.expolist]
    combined_coeffs = [1] # Jacobian coefficient can be ignored when searching for symmetries

    # process `cast`
    for index,prod in enumerate(sector.cast):
        # symmetries may be hidden by the factorization --> undo it
        multiplied_prod = prod.factors[1].copy()
        multiplied_prod.expolist += prod.factors[0].expolist[0]
        multiplied_prod.coeffs *= prod.factors[0].coeffs[0]
        combined_expolists.extend(multiplied_prod.expolist)
        # must distinguish between the individual polynomials and between `Polynomial` and `ExponentiatedPolynomial`
        if type(prod.factors[0]) is ExponentiatedPolynomial:
            combined_coeffs.extend(multiplied_prod.coeffs * sympify_expression('SecDecInternalExponent(%s)*SecDecInternalCast(%i)'%(prod.factors[0].exponent,index)))
        else:
            combined_coeffs.extend(multiplied_prod.coeffs * sympify_expression('SecDecInternalCast(%i)'%index))

    # process `other`
    for index,poly in enumerate(sector.other):
        combined_expolists.append(poly.expolist)
        # must distinguish between `Polynomial` and `ExponentiatedPolynomial`
        if type(poly) is ExponentiatedPolynomial:
            combined_coeffs.extend(poly.coeffs * sympify_expression('SecDecInternalExponent(%s)*SecDecInternalOther(%i)'%(poly.exponent,index)))
        else:
            combined_coeffs.extend(poly.coeffs * sympify_expression('SecDecInternalOther(%i)'%index))

    # return as type `numpy.ndarray`
    return np.vstack(combined_expolists), np.hstack(combined_coeffs)

def _collision_safe_hash(iterable):
    '''
    Return an array containing the hashes of
    an `iterable`. If there is a hash collision
    (if the hashes of two objects that do not
    compare equal are equal), alter the hashes
    until there is no collision any more.

    :param iterable:
        iterable of hashable objects;
        The objects to be hashed.

    '''
    objects = np.asarray(iterable)
    hashes = np.array([hash(i) for i in objects])

    number_of_differing_objects = len( set(objects) )
    number_of_differing_hashes = len( set(hashes) )

    while number_of_differing_objects != number_of_differing_hashes: # while there is a hash collision
        # find the hash with collision
        for i,hash_value in enumerate(hashes):
            # There is a hash collision with ``objects[i]`` if its hash occurs
            # later in the list of `hashes`, but the ``objects[i]`` itself does
            # not appear any more.
            if hash_value in hashes[i+1:] and objects[i] not in objects[i+1:]:
                break

        # add one to the hashes of one of the corresponding object
        hashes[np.where(objects == objects[i])] += 1

        # update `number_of_differing_hashes`
        number_of_differing_hashes = len( set(hashes) )

    return np.array(hashes)

def _array_to_dreadnaut(expolist, coeffs, unique_exponents, unique_coeffs,
                        dreadnaut_file, dreadnaut_hash_file, dreadnaut_canonical_file, workdir):
    '''
    Convert :class:`.Polynomial` `expolist` and `coeffs` to a
    `dreadnaut` graph and write it to a file.

    :param expolist:
        iterable of iterables;
        The variable's powers for each term.

        ..note::
            Each element in the `expolist` will be cast
            to an `int`.

    :param coeffs:
        1d array-like with numerical or sympy-symbolic
        (see http://www.sympy.org/) content, e.g. [x,1,2]
        where x is a sympy symbol;
        The coefficients of the polynomial.

        ..note::
            Each element in the `coeffs` will be cast
            to a `str`.

    :param unique_exponents:
        A 1d array-like with int-like elements;
        An ordered list of all different exponents that
        appear in the input `expolist` or in any `expolist`
        that will be compared to the input `expolist`.

    :param unique_coeffs:
        A 1d array-like with str-like elements;
        An ordered list of all different coefficients that
        appear in the input `coeffs` or in any `coeffs`
        that will be compared to the input `coeffs`.


    :param dreadnaut_file:
        string;
        the filename of the `dreadnaut` input card to be written.

    :param dreadnaut_hash_file:
        string;
        the filename that `dreadnaut` should write the canonical
        graph hash to.

    :param dreadnaut_canonical_file:
        string;
        the filename that `dreadnaut` should write the canonical
        graph to.

    :param workdir:
        string;
        The directory for the communication with `dreadnaut`.
        A directory with the specified name must exist.

        .. note::
           The communication with `dreadnaut` is done via
           files.
    '''
    assert len(expolist.shape) == 2, "_array_to_dreadnaut passed invalid expolist (not 2 x 2)"
    assert expolist.shape[0] > 0, "_array_to_dreadnaut passed invalid expolist (no terms)"
    assert expolist.shape[1] > 0,  "_array_to_dreadnaut passed invalid expolist (no variables)"
    assert len(unique_exponents) > 0, "_array_to_dreadnaut passed no unique exponents"
    assert len(unique_coeffs) > 0, "_array_to_dreadnaut passed no unique coefficients"

    # Note: we can also accept types similar to int such as np.int64
    #    for exponent in unique_exponents:
    #        assert type(exponent) is int, "_array_to_dreadnaut passed non-integer exponents (" + str(type(exponent)) + ")"
    # Note: we can also accept types similar to str such as np.string_
    #    for coeff in unique_coeffs:
    #        assert type(coeff) is str, "_array_to_dreadnaut passed non-string coeffs (" + str(type(coeff)) + ")"

    rows = expolist.shape[0]
    cols = expolist.shape[1]
    elements = rows * cols
    number_unique_exponents = len(unique_exponents)
    number_unique_coeffs = len(unique_coeffs)

    # Number of vertices in graph
    n = cols + rows + elements + number_unique_exponents + number_unique_coeffs

    # Number of vertices generated so far
    offset = 0

    with open(os.path.join(workdir,dreadnaut_file), 'w') as f:

        f.write("n=" + str(n) + "\n")
        f.write("g" + "\n")  # g = begin entering graph

        # Columns should be linked to the elements in steps of cols
        f.write("! columns\n")
        for column in range(0, cols):
            f.write(str(column) + ": " +
                    ",".join(map(str, range(column + cols + rows, column + cols + rows + elements, cols))) + ";" + "\n")
        offset += cols

        # Rows should be linked to the elements in steps of 1
        f.write("! rows\n")
        for row in range(0, rows):
            f.write(str(cols + row) + ": " +
                    ",".join(map(str, range(row * cols + cols + rows, row * cols + rows + cols + cols))) + ";" + "\n")
        offset += rows

        # Create a dictionary with keys given by the unique exponents and entries given by their vertex number
        unique_exponent_count = 0
        exponent_dictionary = dict()
        for exponent in unique_exponents:
            exponent_dictionary[int(exponent)] = offset + elements + unique_exponent_count
            unique_exponent_count += 1

        f.write("! elements\n")
        # Iterate over elements of array row by row
        for x in np.nditer(expolist, order='C'):
            f.write(str(offset) + ": " + str(exponent_dictionary.get(int(x))) + ";" + "\n")
            offset += 1

        f.write("! exponents\n")
        for exponents in range(0, number_unique_exponents):
            f.write(str(offset) + ": ;" + "\n")
            offset += 1

        f.write("! coeffs\n")
        for coeff in unique_coeffs:
            f.write(str(offset) + ": ")
            terms = []
            for row, term_coeff in enumerate(coeffs):
                if str(term_coeff) == coeff:
                    terms.append(cols + row)
            f.write(",".join(map(str, terms)))
            f.write(";" + "\n")
            offset += 1

        # Colour vertices: columns | rows | elements | exponent1 | exponent2 | ... | coeff1 | coeff2 | ...
        f.write("f=[" + \
                str(0) + ":" + str(cols - 1) + "|" + \
                str(cols) + ":" + str(cols + rows - 1) + "|" + \
                str(cols + rows) + ":" + str(cols + rows + elements - 1) + "|" + \
                "|".join(map(str, range(cols + rows + elements,
                                        cols + rows + elements + number_unique_exponents))) + "|" + \
                "|".join(map(str, range(cols + rows + elements + number_unique_exponents,
                                        cols + rows + elements + number_unique_exponents + number_unique_coeffs))) + \
                "]" + "\n")

        # Issue dreadnaut commands
        f.write("c" + "\n")  # c = enable canonical label
        f.write("x" + "\n")  # x = execute
        f.write("B" + "\n")  # B = flush output after every command
        f.write(">" + dreadnaut_hash_file + "\n") # set output to hash file
        f.write("z" + "\n")  # z = print hash (equivalent hashes => graphs possibly the same, inequivalent hashes => graphs not the same)
        f.write(">" + dreadnaut_canonical_file + "\n") # set output to canonical graph file
        f.write("b"+ "\n") # b = print canonical graph (equivalent canonical graphs => equivalent input graphs)
        f.write("->" + "\n") # set output standard out
        f.write("q" + "\n")  # q = quit dreadnaut

def squash_symmetry_redundant_sectors_dreadnaut(sectors, indices=None, dreadnaut='dreadnaut', workdir='dreadnaut_tmp', keep_workdir=False):
    '''
    Reduce a list of sectors by squashing duplicates
    with equal integral.

    Each :class:`.Sector` is converted to a :class:`.Polynomial`
    which is represented as a graph following
    the example of [MP+14]_
    (v2.6 Figure 7, Isotopy of matrices).

    We first multiply each polynomial in the sector by a unique tag
    then sum the polynomials of the sector, this converts a
    sector to a polynomial.
    Next, we convert the `expolist` of the resulting
    polynomial to a graph where each
    unique exponent in the `expolist` is considered to be a
    different symbol.
    Each unique coefficient in the
    polynomial`s `coeffs` is assigned a vertex
    and connected to the row vertex of any term it
    multiplies. The external program `dreadnaut` is then used
    to bring the graph into a canonical form and provide a hash.
    Sectors with equivalent hashes may be identical,
    their canonical graphs are compared and if they are identical
    the sectors are combined.

    .. note::
        This function calls the command line executable of
        `dreadnaut` [MP+14]_.
        It has been tested with `dreadnaut` version nauty26r7.

    See also:
    :func:`squash_symmetry_redundant_sectors_sort`

    :param sectors:
        iterable of :class:`.Sector`; the sectors to be
        reduced.

    :param indices:
        iterable of integers, optional;
        The indices of the variables to consider. If not
        provided, all indices are taken into account.

    :param dreadnaut:
        string;
        The shell command to run `dreadnaut`.

    :param workdir:
        string;
        The directory for the communication with `dreadnaut`.
        A directory with the specified name will be created
        in the current working directory. If the specified
        directory name already exists, an :class:`OSError`
        is raised.

        .. note::
            The communication with `dreadnaut` is done via
            files.

    :param keep_workdir:
        bool;
        Whether or not to delete the `workdir` after execution.

    '''
    if not isinstance(sectors, list):
        sectors = list(sectors)

    if indices is None:
        replaced_sectors = sectors
    else:
        # must move the indices to ignore to the coefficients
        replaced_sectors = _remove_variables(sectors, indices)

    os.mkdir(workdir)
    try:
        # create list of all exponents that appear in any polynomial in any sector
        # create list of all coefficients that appear in any polynomial in any sector
        unique_exponents = []
        unique_coeffs = []
        for sector in replaced_sectors:
            this_sector_expolist, this_sector_coeffs = _sector2array(sector)

            # Note:
            # - we can not use _collision_safe_hash as the hashes may collide between sectors
            # - we must use str as np.unique does not function correctly when array elements can not be sorted (e.g. sympy symbols)
            this_sector_coeffs = [str(coeff) for coeff in this_sector_coeffs]

            unique_exponents += np.unique(this_sector_expolist).tolist()
            unique_coeffs += np.unique(this_sector_coeffs).tolist()

        unique_exponents = np.unique(unique_exponents)
        unique_coeffs = np.unique(unique_coeffs)

        dreadnaut_hash_file_suffix = "_hash"
        dreadnaut_canonical_file_suffix = "_canonical"

        sector_number = 0
        sector_hash_array = []
        for sector in replaced_sectors:
            this_sector_expolist, this_sector_coeffs = _sector2array(sector)

            sector_name = str(sector_number)
            dreadnaut_file = sector_name
            dreadnaut_hash_file = sector_name + dreadnaut_hash_file_suffix
            dreadnaut_canonical_file = sector_name + dreadnaut_canonical_file_suffix
            command_line_command = dreadnaut
            _array_to_dreadnaut(this_sector_expolist, this_sector_coeffs, unique_exponents, unique_coeffs,\
                                dreadnaut_file, dreadnaut_hash_file, dreadnaut_canonical_file, workdir)

            # redirect dreadnaut stdout
            with open(os.path.join(workdir, 'stdout'), 'w') as stdout:
                # redirect dreadnaut stderr
                with open(os.path.join(workdir, 'stderr'), 'w') as stderr:
                    # redirect dreadnaut stdin
                    with open(os.path.join(workdir, dreadnaut_file), 'r') as stdin:
                        # run dreadnaut
                        #    subprocess.check_call --> run dreadnaut, block until it finishes and raise error on nonzero exit status
                        try:
                            subprocess.check_call(command_line_command, stdin=stdin,
                                                  stdout=stdout, stderr=stderr, cwd=workdir)
                        except OSError as error:
                            if dreadnaut not in str(error):
                                error.filename = dreadnaut
                            raise

            # collect dreadnaut canonical hash for sector
            with open(os.path.join(workdir, dreadnaut_hash_file), 'r') as f:
                sector_hash = f.readline()
            assert len(sector_hash) > 0, "dreadnaut hash file (" + dreadnaut_hash_file + ") empty"

            sector_hash_array.append(sector_hash)
            sector_number += 1

        sector_sort_index = np.argsort(sector_hash_array)

        # Note:
        # - if sector hashes are not equal the sectors are different
        # - if sector hashes are equal then the sectors could be equivalent but this is not guaranteed,
        #   we must check that the canonical graphs are identical

        # iterate through the list of sectors, if the hashes are the same check that the canonical graph is the same, if so, squash
        # since we know the sorted indices of the sector's hashes we only need to consider the sectors with consecutive indices
        previous_index = sector_sort_index[0]
        previous_sector_hash = sector_hash_array[previous_index]
        previous_sector = sectors[previous_index].copy()
        output = [previous_sector]
        for sector_index in sector_sort_index[1:]:
            if sector_hash_array[sector_index] == previous_sector_hash:
                # load the canonical graphs are really equal
                with open(os.path.join(workdir, str(previous_index)+dreadnaut_canonical_file_suffix)) as f:
                    previous_canonical_graph = f.readlines()
                with open(os.path.join(workdir, str(sector_index)+dreadnaut_canonical_file_suffix)) as f:
                    canonical_graph = f.readlines()

                # The first few lines of the canonical graph file give the permutation of the vertices
                # required to reach canonical form, this can differ for isomorphic graphs so we ignore it
                # search for the first vertex which has a line beginning "  0 : "
                first_line_previous_canonical_graph = 0
                for line in previous_canonical_graph:
                    if line[:6] == "  0 : ":
                        break
                    first_line_previous_canonical_graph += 1
                first_line_canonical_graph = 0
                for line in canonical_graph:
                    if line[:6] == "  0 : ":
                        break
                    first_line_canonical_graph += 1

                if canonical_graph[first_line_canonical_graph:] == \
                        previous_canonical_graph[first_line_previous_canonical_graph:]:
                    # sector is equal to the previous sector --> squash into previous sector
                    sector_is_new = False
                    # squash this sector into the previous one
                    previous_sector.Jacobian.coeffs[0] += sectors[sector_index].Jacobian.coeffs[0]
                else:
                    # hash collision, this sector is not equal to the previous one --> update output
                    sector_is_new = True
            else:
                # hash does not match previous sector, sector is a new sector
                sector_is_new = True

            if sector_is_new:
                # this sector is not equal to the previous one --> update output
                previous_index = sector_index
                previous_sector_hash = sector_hash_array[sector_index]
                previous_sector = sectors[sector_index].copy()
                output.append(previous_sector)

    finally:
        if not keep_workdir:
            shutil.rmtree(workdir)

    return output

def _remove_variables(sectors, indices_to_keep):
    '''
    Remove the indices of the variables absent in
    `indices_to_keep` from `sectors`.

    '''
    first_sector = sectors[0]
    all_indices = set(range(first_sector.number_of_variables))
    indices_to_remove = sorted(all_indices - set(indices_to_keep), reverse=True)
    symbols = first_sector.Jacobian.symbols

    def remove_indices(expression):
        new_expression = expression
        for index in indices_to_remove:
            new_expression = new_expression.replace(index, symbols[index], True) # last argument is `remove`
        return new_expression

    output = []
    for sector in sectors:
        new_Jacobian = remove_indices(sector.Jacobian)
        new_cast = [remove_indices(prod) for prod in sector.cast]
        new_other = [remove_indices(poly) for poly in sector.other]

        output.append( Sector(new_cast, new_other, new_Jacobian) )

    return output

def squash_symmetry_redundant_sectors_sort(sectors, sort_function, indices=None):
    '''
    Reduce a list of sectors by squashing duplicates
    with equal integral.

    If two sectors only differ by a permutation of the
    polysymbols (to be interpreted as integration
    variables over some inteval), then the two sectors
    integrate to the same value. Thus we can drop one
    of them and count the other twice. The multiple
    counting of a sector is accounted for by increasing
    the coefficient of the Jacobian by one.

    Equivalence up to permutation is established by
    applying the `sort_function` to each sector,
    this brings them into a canonical form.
    Sectors with identical canonical forms differ only
    by a permutation.

    Note: whether all symmetries are found depends on
    the choice of `sort_function`. The sort function
    :func:`pySecDec.matrix_sort.Pak_sort` should find
    all symmetries whilst the sort functions
    :func:`pySecDec.matrix_sort.iterative_sort`
    and
    :func:`pySecDec.matrix_sort.light_Pak_sort` are
    faster but do not identify all symmetries.

    See also:
    :func:`squash_symmetry_redundant_sectors_dreadnaut`

    Example:

    >>> from pySecDec.algebra import Polynomial
    >>> from pySecDec.decomposition import Sector
    >>> from pySecDec.decomposition import squash_symmetry_redundant_sectors_sort
    >>> from pySecDec.matrix_sort import Pak_sort
    >>>
    >>> poly = Polynomial([(0,1),(1,0)], ['a','b'])
    >>> swap = Polynomial([(1,0),(0,1)], ['a','b'])
    >>> Jacobian_poly = Polynomial([(1,0)], [3]) # three
    >>> Jacobian_swap = Polynomial([(0,1)], [5]) # five
    >>> sectors = (
    ...               Sector([poly],Jacobian=Jacobian_poly),
    ...               Sector([swap],Jacobian=Jacobian_swap)
    ...           )
    >>>
    >>> reduced_sectors = squash_symmetry_redundant_sectors_sort(sectors,
    ...                   Pak_sort)
    >>> len(reduced_sectors) # symmetry x0 <--> x1
    1
    >>> # The Jacobians are added together to account
    >>> # for the double counting of the sector.
    >>> reduced_sectors[0].Jacobian
     + (8)*x0

    :param sectors:
        iterable of :class:`.Sector`; the sectors to be
        reduced.

    :param sort_function:
        :func:`pySecDec.matrix_sort.iterative_sort`,
        :func:`pySecDec.matrix_sort.light_Pak_sort`, or
        :func:`pySecDec.matrix_sort.Pak_sort`;
        The function to be used for finding a canonical
        form of the sectors.

    :param indices:
        iterable of integers, optional;
        The indices of the variables to consider. If not
        provided, all indices are taken into account.

    '''
    if not isinstance(sectors, list):
        sectors = list(sectors)

    if indices is None:
        replaced_sectors = sectors
    else:
        # must move the indices to ignore to the coefficients
        replaced_sectors = _remove_variables(sectors, indices)

    # combine all expolists and coeffs into two large arrays
    all_sectors_expolist = []
    all_sectors_coeffs = []
    for sector in replaced_sectors:
        this_sector_expolist, this_sector_coeffs = _sector2array(sector)

        all_sectors_expolist.append( this_sector_expolist )
        all_sectors_coeffs.append( this_sector_coeffs )

    # call the collision free hash for the coefficients
    # must call the collision safe hash on ALL coefficients TOGETHER
    all_sectors_coeffs = np.array(all_sectors_coeffs)
    all_sectors_coeffs = _collision_safe_hash(all_sectors_coeffs.flatten()).reshape(all_sectors_coeffs.shape)

    # combine coefficients and expolists into one large array
    assert len(all_sectors_expolist)  == len(all_sectors_coeffs)
    all_sectors_array = []
    for this_sector_expolist,this_sector_coeffs in zip(all_sectors_expolist,all_sectors_coeffs):
        this_sector_array = np.hstack((this_sector_coeffs.reshape(-1,1),this_sector_expolist))

        # sort the arrays of the individual sectors to pick one specific permutation --> symmetry finding
        sort_function(this_sector_array)

        all_sectors_array.append( this_sector_array )

    # clean up large temporary arrays
    del all_sectors_expolist, all_sectors_coeffs

    all_sectors_array = np.array(all_sectors_array)

    # sort the sectors
    indices_sorted_sectors = argsort_ND_array(all_sectors_array)

    # prune sectors that are equal
    # Since we sorted the sectors above, we only need to consider the
    # sectors that next to each other.
    previous_index = indices_sorted_sectors[0]
    previous_sector_array = all_sectors_array[previous_index]
    previous_sector = sectors[previous_index].copy()
    output = [previous_sector]
    for sector_index in indices_sorted_sectors[1:]:
        if np.array_equal(all_sectors_array[sector_index],previous_sector_array):
            # squash this sector into the previous one
            previous_sector.Jacobian.coeffs[0] += sectors[sector_index].Jacobian.coeffs[0]
        else:
            # this sector is not equal to the previous one --> update output
            previous_index = indices_sorted_sectors[sector_index]
            previous_sector_array = all_sectors_array[sector_index]
            previous_sector = sectors[sector_index].copy()
            output.append(previous_sector)

    return output
