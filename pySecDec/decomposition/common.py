'''
This module collects routines that are used by
multiple decompition modules.

'''

from ..algebra import Polynomial, ExponentiatedPolynomial, Product
from ..matrix_sort import iterative_sort, Pak_sort
from ..misc import argsort_ND_array
import numpy as np
import sympy as sp

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

        if Jacobian is not None:
            assert len(Jacobian.coeffs) == 1, "`Jacobian` must be a monomial"
        else:
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

        self.cast = []
        for item in cast:
            if hasattr(item, 'factors'): # expect type `Product`
                self.cast.append(item.copy())
            else: # expect type `Polynomial`
                if hasattr(item, 'exponent'): # `ExponentiatedPolynomial` --> same exponent for the factorizable monomial
                    monomial = ExponentiatedPolynomial(initial_monomial_factor.expolist,
                                                       initial_monomial_factor.coeffs,
                                                       polysymbols=initial_monomial_factor.polysymbols,
                                                       exponent=item.exponent)
                    item = Product(monomial, item)
                else:
                    item = Product(initial_monomial_factor, item)

                # try to factorize only if constructed from `Polynomial`
                refactorize(item)
                self.cast.append(item)

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
        factors_as_polynomials = [Polynomial(factor.expolist, factor.coeffs, factor.polysymbols, copy=False) for factor in prod.factors]
        multiplied_prod = factors_as_polynomials[0] * factors_as_polynomials[1]
        combined_expolists.extend(multiplied_prod.expolist)
        # must distinguish between the individual polynomials and between `Polynomial` and `ExponentiatedPolynomial`
        if type(prod.factors[0]) is ExponentiatedPolynomial:
            combined_coeffs.extend(multiplied_prod.coeffs * sp.sympify('SecDecInternalExponent(%s)*SecDecInternalCast(%i)'%(prod.factors[0].exponent,index)))
        else:
            combined_coeffs.extend(multiplied_prod.coeffs * sp.sympify('SecDecInternalCast(%i)'%index))

    # process `other`
    for index,poly in enumerate(sector.other):
        combined_expolists.append(poly.expolist)
        # must distinguish between `Polynomial` and `ExponentiatedPolynomial`
        if type(prod.factors[0]) is ExponentiatedPolynomial:
            combined_coeffs.extend(poly.coeffs * sp.sympify('SecDecInternalExponent(%s)*SecDecInternalOther(%i)'%(prod.factors[0].exponent,index)))
        else:
            combined_coeffs.extend(poly.coeffs * sp.sympify('SecDecInternalOther(%i)'%index))

    # return as type `numpy.ndarray`
    return np.vstack(combined_expolists), np.hstack(combined_coeffs)

def _collision_safe_hash(iterable):
    '''
    Return an array containingh the hashes of
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

def squash_symmetry_redundant_sectors(sectors, sort_function):
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

    Example:

    >>> from pySecDec.algebra import Polynomial
    >>> from pySecDec.decomposition import Sector
    >>> from pySecDec.decomposition import squash_symmetry_redundant_sectors
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
    >>> reduced_sectors = squash_symmetry_redundant_sectors(sectors,
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
        :func:`pySecDec.matrix_sort.iterative_sort` or
        :func:`pySecDec.matrix_sort.Pak_sort`;
        The function to be used for finding a canonical
        form of the sectors.

    '''
    if not isinstance(sectors, list):
        sectors = list(sectors)

    # combine all expolists and coeffs into one large array
    # use the collision free hash for the coefficients
    all_sectors_array = []
    for sector in sectors:
        this_sector_expolist, this_sector_coeffs = _sector2array(sector)
        this_sector_coeffs = _collision_safe_hash(this_sector_coeffs)

        this_sector_array = np.hstack([this_sector_coeffs.reshape(len(this_sector_coeffs),1),this_sector_expolist])

        # sort the arrays of the individual sectors to pick one specific permutation --> symmetry finding
        sort_function(this_sector_array)

        all_sectors_array.append( this_sector_array )

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
