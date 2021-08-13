"""

The iterative sector decomposition routines.

"""

from .common import Sector
from ..algebra import Product, refactorize
from ..misc import powerset
import numpy as np

# ********************** primary decomposition **********************

def primary_decomposition_polynomial(polynomial, indices=None):
    r'''
    Perform the primary decomposition on a single polynomial.

    .. seealso::
        :func:`.primary_decomposition`

    :param polynomial:
        :class:`.algebra.Polynomial`;
        The polynomial to eliminate the Dirac delta from.

    :param indices:
        iterable of integers or None;
        The indices of the parameters to be considered as
        integration variables. By default (``indices=None``),
        all parameters are considered as integration
        variables.

    '''
    primary_sectors = []

    # consider all indices if `indices` is None
    if indices is None:
        indices = range(polynomial.expolist.shape[1])

    coeffs = polynomial.coeffs
    expolist = polynomial.expolist
    polysymbols = polynomial.polysymbols

    if len(indices) == 1:
        # do not remove the only Feynman parameter from the `polysymbols`
        sectorpoly = polynomial.copy()
        sectorpoly.expolist[:,indices[0]] = 0

        return [sectorpoly]

    else:
        for i in indices:
            # keep the type (`polynomial` can have a subtype of `Polynomial`)
            sectorpoly = polynomial.copy()

            # "pinch" (delete) Feynman parameter `i`
            #   => this is equivalent to setting the exponent to zero
            #     => that is however equivalent to setting the parameter to one
            sectorpoly.expolist = np.delete(expolist,i,axis=1)
            sectorpoly.polysymbols = np.delete(polysymbols,i)
            sectorpoly.number_of_variables -= 1

            primary_sectors.append(sectorpoly)

    return primary_sectors

def primary_decomposition(sector, indices=None):
    r'''
    Perform the primary decomposition as described in
    chapter 3.2 (part I) of arXiv:0803.4177v2 [Hei08]_.
    Return a list of :class:`.Sector` - the primary
    sectors.
    For `N` Feynman parameters, there are `N` primary
    sectors where the `i`-th Feynman parameter is set to
    `1` in sector `i`.

    .. seealso::
        :func:`.primary_decomposition_polynomial`

    :param sector:
        :class:`.Sector`;
        The container holding the polynomials (typically
        :math:`U` and :math:`F`) to eliminate the Dirac
        delta from.

    :param indices:
        iterable of integers or None;
        The indices of the parameters to be considered as
        integration variables. By default (``indices=None``),
        all parameters are considered as integration
        variables.

    '''
    # consider all variables if `indices` is None
    # convert `indices` to list otherwise
    if indices is None:
        N = sector.number_of_variables
        indices = range(N)
    else:
        indices = list(indices)
        N = len(indices)

    # Must perform `primary_decomposition_polynomial` for each element
    # of `sector.other`, `[sector.Jacobian]`, and each of the two factors
    # of `sector.cast`.
    primary_sector_polynomials_other = [primary_decomposition_polynomial(poly, indices) for poly in sector.other]
    primary_sector_polynomials_Jacobian = primary_decomposition_polynomial(sector.Jacobian, indices)
    primary_sector_polynomials_cast_factor0 = [primary_decomposition_polynomial(polyprod.factors[0], indices) for polyprod in sector.cast]
    primary_sector_polynomials_cast_factor1 = [primary_decomposition_polynomial(polyprod.factors[1], indices) for polyprod in sector.cast]

    # Collect the primary decomposed polynomials back into the `Sector` container class
    primary_sectors = [
        Sector(
                [Product(cast0[sector_index],cast1[sector_index]) for cast0,cast1 in zip(primary_sector_polynomials_cast_factor0, primary_sector_polynomials_cast_factor1)],
                [other[sector_index] for other in primary_sector_polynomials_other],
                primary_sector_polynomials_Jacobian[sector_index]
        ) for sector_index in range(N)
    ]
    return primary_sectors

# ********************** iterative decomposition **********************

class EndOfDecomposition(Exception):
    '''
    This exception is raised if the function
    :func:`.iteration_step` is called although
    the sector is already in standard form.

    '''
    pass

def remap_parameters(singular_parameters, Jacobian, *polynomials):
    r'''
    Remap the Feynman parameters according to eq. (16) of
    arXiv:0803.4177v2 [Hei08]_. The parameter whose index comes
    first in `singular_parameters` is kept fix.

    The remapping is done in place; i.e. the `polynomials` are
    **NOT** copied.

    :param singular_parameters:
        list of integers;
        The indices :math:`\alpha_r` such that at least one
        of `polynomials` becomes zero if all
        :math:`t_{\alpha_r} \rightarrow 0`.

    :param Jacobian:
        :class:`.Polynomial`;
        The Jacobian determinant is multiplied to this polynomial.

    :param polynomials:
        abritrarily many instances of :class:`.algebra.Polynomial`
        where all of these have an equal number of variables;
        The polynomials of Feynman parameters to be remapped. These
        are typically :math:`F` and :math:`U`.

    Example:

    .. code-block:: python

        remap_parameters([1,2], Jacobian, F, U)

    '''
    assert polynomials, "No polynomial for modification passed"

    num_parameters = polynomials[0].expolist.shape[1]

    for poly in polynomials:
        assert num_parameters == poly.expolist.shape[1], 'All polynomials must have the same number of variables'

    for poly in polynomials:
        for param in singular_parameters[1:]:
            poly.expolist[:,singular_parameters[0]] += poly.expolist[:,param] # This modifies in place!

    Jacobian.expolist[:,singular_parameters[0]] += len(singular_parameters) - 1

def find_singular_set(sector, indices=None):
    '''
    Function within the iterative sector decomposition procedure
    which heuristically chooses an optimal decomposition set.
    The strategy was introduced in arXiv:hep-ph/0004013 [BH00]_
    and is described in 4.2.2 of arXiv:1410.7939 [Bor14]_.
    Return a list of indices.

    :param sector:
        :class:`.Sector`;
        The sector to be decomposed.

    :param indices:
        iterable of integers or None;
        The indices of the parameters to be considered as
        integration variables. By default (``indices=None``),
        all parameters are considered as integration
        variables.

    '''
    # consider all indices if `indices` is None
    if indices is None:
        indices = range(sector.number_of_variables)

    def get_poly_to_transform(sector, indices):
        '''
        Return a :class:`.algebra.Product` in
        `sector.cast` that is not in the desired form
        `<monomial> * <const + ...>` yet.
        Raise `EndOfDecomposition` if the desired form is
        reached.

        '''
        for polyprod in sector.cast:
            if not polyprod.factors[1].has_constant_term(indices):
                return polyprod
        # Control only reaches this point if the desired form is
        # already reached for all polynomials in ``sector.cast``.
        raise EndOfDecomposition()

    # find a polynomial to cast that is not in the form ``const + ... yet``
    polyprod = get_poly_to_transform(sector, indices)
    poly = polyprod.factors[1]
    possible_sets = []

    # find sets that nullyfy the selected polynomial
    # only consider sets of the smallest possible size
    for singular_set in powerset(indices,min_length=2):
        if possible_sets and len(possible_sets[0])<len(singular_set):
            break
        if poly.becomes_zero_for(singular_set):
            possible_sets.append(singular_set)
    assert possible_sets

    # First check how many poynomials of the sector nullify for
    # each set of the `possible sets` of fixed length.
    # Second only gather those sets which nullify the most polynomials
    howmany_max = 1
    howmany = 0
    best_sets = []
    for singular_set in possible_sets:
        for polyprod in sector.cast:
            if not polyprod.factors[1].has_constant_term(indices):
                howmany += 1
            if howmany > howmany_max:
                best_sets = []
                howmany_max = howmany
        best_sets.append(singular_set)
        howmany = 0

    # Choose the set of Feynman parameters with the
    # highest powers for remapping.
    exposum_max = -np.inf
    for test_set in best_sets:
        exposum = poly.expolist[:,test_set].max(axis=0).sum()
        if exposum > exposum_max:
            exposum_max = exposum
            best_set = test_set
    assert np.isfinite(exposum_max)

    # return the chosen set
    return best_set

def iteration_step(sector, indices=None):
    '''
    Run a single step of the iterative sector decomposition as described
    in chapter 3.2 (part II) of arXiv:0803.4177v2 [Hei08]_.
    Return an iterator of :class:`.Sector` - the arising subsectors.

    :param sector:
        :class:`.Sector`;
        The sector to be decomposed.

    :param indices:
        iterable of integers or None;
        The indices of the parameters to be considered as
        integration variables. By default (``indices=None``),
        all parameters are considered as integration
        variables.

    '''
    # consider all indices if `indices` is None
    if indices is None:
        indices = range(sector.number_of_variables)

    # find a set that describes the transformation to be performed
    singular_set = find_singular_set(sector, indices)

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
        polynomials_to_transform = [polyprod.factors[0] for polyprod in subsector.cast] + \
                                   [polyprod.factors[1] for polyprod in subsector.cast] + \
                                   subsector.other + [subsector.Jacobian]
        remap_parameters(singular_set, subsector.Jacobian, *polynomials_to_transform)
        for polyprod in subsector.cast:
            refactorize(polyprod,singular_set[0])
        yield subsector

def iterative_decomposition(sector, indices=None):
    '''
    Run the iterative sector decomposition as described
    in chapter 3.2 (part II) of arXiv:0803.4177v2 [Hei08]_.
    Return an iterator of :class:`.Sector` - the
    arising subsectors.

    :param sector:
        :class:`.Sector`;
        The sector to be decomposed.

    :param indices:
        iterable of integers or None;
        The indices of the parameters to be considered as
        integration variables. By default (``indices=None``),
        all parameters are considered as integration
        variables.

    '''
    # convert `indices` to list if not None
    if indices is not None:
        indices = list(indices)

    try:
        subsectors = iteration_step(sector, indices) # only this line can raise `EndOfDecomposition`
        for subsector in subsectors:
            for deeper_subsector in iterative_decomposition(subsector, indices):
                yield deeper_subsector
    except EndOfDecomposition:
        yield sector
