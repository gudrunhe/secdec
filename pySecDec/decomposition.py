"""The sector decomposition routines"""

from .polynomial import Polynomial
import numpy as np

def decompose_step(singular_parameters, Jacobian, *polynomials):
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

        decompose_step([1,2], Jacobian, F, U)

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
