"""The sector decomposition routines"""

from .polynomial import Polynomial, PolynomialProduct
from .sector import Sector
import numpy as np

# ********************** primary decomposition **********************

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

# ********************** iterative decomposition **********************

class EndOfDecomposition(Exception):
    '''
    This exception is raised if the function
    :func:`.iterative_step` is called although
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
    assert Jacobian.coeffs == [""], "`Jacobian` must be a monomial without coefficient"
    assert polynomials, "No polynomial for modification passed"

    num_parameters = polynomials[0].expolist.shape[1]

    for poly in polynomials:
        assert num_parameters == poly.expolist.shape[1], 'All polynomials must have the same number of variables'

    for poly in polynomials:
        for param in singular_parameters[1:]:
            poly.expolist[:,singular_parameters[0]] += poly.expolist[:,param] # This modifies in place!

    Jacobian.expolist[:,singular_parameters[0]] += len(singular_parameters) - 1
