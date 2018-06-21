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

def apply_regions( exp_param , sector, normaliz='normaliz', workdir='normaliz_tmp'):
    '''
    Unexpanded integrals for expansion by regions in limit `exp_param` to zero.

    .. note::
        This function calls the command line executable of
        `normaliz` [BIR]_. See :ref:`installation_normaliz`
        for installation and a list of tested versions.

    :param exp_param:
        string;
        Expansion parameter.

    :param sector:
        :class:`.Sector`;
        The sector to be decomposed into regions.

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
    # include expansion parameter as polynomial variable
    symbols = [ exp_param ] + sector.Jacobian.polysymbols
    sector.Jacobian = Polynomial.from_expression( str( sector.Jacobian ), symbols )
    for product in sector.cast:
        for i, factor in enumerate( product.factors ):
            product.factors[i] = Polynomial.from_expression( str( factor ), symbols )
    for i, poly in enumerate( sector.other ):
        sector.other[i] = Polynomial.from_expression( str( poly ), symbols )

    # should sector.other be included in the decomposition?
    region_vectors = find_regions( symbols.index( exp_param ), [ prod.factors[1] for prod in sector.cast ] + [ prod.factors[1] for prod in sector.other ] , normaliz, workdir )

    def rescale_polynomial(polynomial, region_vector):
        return Polynomial([np.append( exponentvector,  np.dot(region_vector, exponentvector)) for exponentvector in polynomial.expolist ], polynomial.coeffs, polysymbols = symbols + ['zz'])

    def make_sector_region(region_vector):
        subsector = sector.copy()
        # includes scaling factor from integral measure
        subsector.Jacobian =  Polynomial( [np.append( np.zeros(len(symbols), dtype = int), sum(region_vector[1:]))], [1], polysymbols = symbols + ['zz'] )*rescale_polynomial( subsector.Jacobian, region_vector )
        for resulting_product, output_product in zip(sector.cast, subsector.cast):
            for j in range(2):
                output_product.factors[j] = rescale_polynomial( resulting_product.factors[j], region_vector )
            refactorize(output_product)
        for resulting_polynomial, output_polynomial in zip(sector.other, subsector.other):
            output_polynomial = rescale_polynomial(resulting_polynomial)
        return subsector

    #ToDo: pass relation between 'zz' and exp_param (i.e. exp_param = zz**somepower)
    for region in region_vectors:
        yield make_sector_region( region )
