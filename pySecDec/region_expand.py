"""
Expansion by Regions
--------------------

Routines to perform an expansion by regions, see e.g. [PS11]_.

"""

from .decomposition.common import Sector, refactorize
from .polytope import Polytope
from .algebra import Polynomial, ExponentiatedPolynomial
import numpy as np, sympy as sp

def find_regions( exp_param_index , polynomials, normaliz='normaliz', workdir='normaliz_tmp'):
    '''
    Find regions for the expansion by regions
    as described in [PS11]_.

    .. note::
        This function calls the command line executable of
        `normaliz` [BIR]_. See :ref:`installation_normaliz`
        for installation and a list of tested versions.

    :param exp_param_index:
        int;
        The index of the expansion parameter in the expolist.

    :param polynomials:
        list of instances of :class:`.Polynomial` where
        all of these have an equal number of variables;
        The polynomials to calculate the regions for.

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
    sum = 0
    for poly in polynomials:
        sum += Polynomial(poly.expolist, np.ones_like(poly.coeffs, dtype=int))

    polytope_vertices = sum.expolist
    polytope = Polytope(vertices=polytope_vertices)
    polytope.complete_representation(normaliz, workdir)

    facets = polytope.facets[:,:-1] # do not need offset term "a_F"

    regions = facets[ facets[:,exp_param_index] > 0 ]

    return regions

def apply_region(polynomials, region_vector, expansion_parameter_index):
    r'''
    Apply the `region_vector` to the input `polynomials`.

    .. note::
        Redefines the `expansion_parameter` as
        :math:`\rho \rightarrow \rho ^ n`, where
        :math:`n` is given by the `region_vector`.

    .. note::
        `apply_region` modifies the input
        `polynomials`.

    :param polynomials:
        iterable of polynomials;
        Polynomials to be computed in different regions.

    :param region_vector:
        vector-like array;
        Region specified by the power of the `expansion_parameter`
        in the rescaled variables. The region vectors have to be
        specified in the same order as the symbols are specified
        in the polynomials. E.g. if symbols are specified as
        ['x0','x1','rho'] and want rescaling x0 -> rho^i * x0,
        x1 -> rho^k * x1 and rho -> rho^n, then the region vector
        needs to be [i,k,n]

    :param expansion_parameter_index:
        integer;
        Index of the expansion parameter in the list of symbols.

    '''
    for polynomial in polynomials:
        polynomial.expolist[:,expansion_parameter_index] = np.dot(polynomial.expolist, region_vector)

    return polynomials


def derive_prod(poly_list,numerator,index,polynomial_name_indices):
    r"""
    Calculates the derivative of a product of polynomials using

    .. math::
        \frac{\partial}{\partial x_i} \left( \prod_j P_j^{\alpha_j} N \right)
        = \prod_j P_j^{\alpha_j-1} N'

    where N' is given by

    .. math::
        N' = \left( \sum_j N \alpha_j \frac{\partial P_j}{\partial x_i} \prod_{k \neq j} P_k \right)
        + \left(\prod_l P_l \right) \left[ \left( \sum_k \frac{\partial N}{\partial P_k}
        \frac{\partial P_k}{\partial x_i} \right) + \frac{\partial N}{\partial x_i} \right] .

    :param poly_list:
        list of :class:`.ExponentiatedPolynomial`;
        The exponentiated polynomials that should be differentiated.
        They need to be defined in terms of the symbols ``x0,x1,x2,..`` and
        ``p0,p1,p2..`` where ``p0,p1,p2..`` are the bases of the exponentiated
        polynomials.

    :param numerator:
        :class:`.Polynomial`;
        The numerator also defined as an exponentiated polynomial with
        ``symbols = [x0,x1,...,p0,p1,...]``.

    :param index:
        integer;
        Index of variable with respect to which the derivative is taken.

    :param polynomial_name_indices:
        iterable;
        Indices of polynomial names in poly_symbols.


    """
    number_of_polys = len(poly_list)
    symbols = numerator.polysymbols
    number_of_symbols = len(numerator.polysymbols)

    explicit_poly_list = []
    dummy_poly_list = []
    for i in range(number_of_polys):
        # define the polynomials without exponents because needed for explicit derivative
        explicit_poly = Polynomial(poly_list[i].expolist, poly_list[i].coeffs,poly_list[i].polysymbols)
        explicit_poly_list.append(explicit_poly)

        # define the polynomials as just 1*p1, 1*p2 etc. as only want to evaluate derivatives explicitly
        expolist = np.zeros(number_of_symbols,dtype=np.int)
        expolist[polynomial_name_indices[i]] = 1
        dummy_poly = Polynomial([expolist], [1], symbols)
        dummy_poly_list.append(dummy_poly)

    # calculate the product of all polynomials
    prod_all_polys = np.prod(dummy_poly_list)

    # calculate the new numerator according to formula:
    numerator_new = numerator.derive(index)*prod_all_polys

    for i in range(number_of_polys):
        prod = 1
        for j in range(number_of_polys):

            if j != i:

                prod *= dummy_poly_list[j]

        summand1 = numerator * poly_list[i].exponent * explicit_poly_list[i].derive(index) * prod
        summand2 = numerator.derive(polynomial_name_indices[i]) * explicit_poly_list[i].derive(index) * prod_all_polys

        numerator_new += (summand1+summand2)

    new_poly_list=[]
    for i in range(number_of_polys):

        # TODO: what to do if power becomes 0 ??
        new_poly_list.append(ExponentiatedPolynomial(poly_list[i].expolist, poly_list[i].coeffs,poly_list[i].exponent -1 ,poly_list[i].polysymbols))

    return new_poly_list,numerator_new
