"""
This module implements a wrapper for loop integrals
around the main program
:func:`pySecDec.code_writer.make_package` - the
function :func:`.loop_package`.

"""

from ..code_writer import make_package
import numpy as np
import sympy as sp

def loop_package(target_directory, name, loop_integral, requested_order,
                 real_parameters=[], complex_parameters=[],
                 contour_deformation=True, additional_prefactor=1,
                 form_optimization_level=2, form_work_space='500M',
                 stabilize=False, decomposition_method='iterative'):
    '''
    Decompose, subtract and expand a Feynman
    parametrized loop integral. Return it as
    c++ package.

    .. note::
        When instanciating the
        :class:`pySecDec.loop_integral.LoopIntegral`,
        use the symbol ``I`` for the imaginary unit.

    .. seealso::
        This function is a wrapper around
        :func:`pySecDec.code_writer.make_package`.

    :param target_directory:
        string;
        The output directory.

    :param name:
        string;
        The name of the c++ namespace.

    :param loop_integral:
        :class:`pySecDec.loop_integral.LoopIntegral`;
        The loop integral to be computed.

    :param requested_orders:
        integer;
        Compute the expansion in the regulator to this
        order.

    :param real_parameters:
        iterable of strings or sympy symbols, optional;
        Parameters to be interpreted as real numbers,
        e.g. Mandelstam invariants and masses.

    :param complex parameters:
        iterable of strings or sympy symbols, optional;
        Parameters to be interpreted as complex numbers.
        To use the complex mass scheme, define the masses
        as complex parameters.

    :param contour_deformation:
        bool, optional;
        Whether or not to produce code for contour
        deformation.
        Default: ``True``.

    :param additional_prefactor:
        string or sympy expression;
        An additional factor to be multiplied to the loop
        integral. It may depend on the regulator, the
        `real_parameters`, and the `complex_parameters`.

    :param form_optimization_level:
        integer out of the interval [0,3], optional;
        The optimization level to be used in FORM.
        Default: ``2``.

    :param form_work_space:
        string, optional;
        The FORM WorkSpace. Default: ``'500M'``.

    :param stabilize:
        bool, optional;
        Whether or not to bring subexpression to a common
        denominator. Default: ``False``.

        .. warning::
            This is very extensive concerning both - the
            algebra and the numerics. It should only be
            set to ``True`` if numerical instabilities
            occur.

    :param decomposition_method:
        string, optional;
        The strategy to decompose the polynomials. The
        following strategies are available:

        * 'iterative' (default)
        * 'geometric'

        .. note::
            For 'geometric', the third-party program "normaliz"
            is needed. See :ref:`installation_normaliz`.

    '''
    # append the regulator to the `polynomials_to_decompose` (`F` and `U`)
    polynomials_to_decompose = [loop_integral.exponentiated_F.copy(), loop_integral.exponentiated_U.copy()]
    for poly in polynomials_to_decompose:
        poly.polysymbols.append(loop_integral.regulator)
        poly.expolist = np.hstack([poly.expolist, np.zeros([len(poly.expolist),1], dtype=int)])

    other_polynomials = [loop_integral.numerator]
    if sp.sympify( loop_integral.measure ) != 1:
        # need ``loop_integral.measure`` only if it is nontrivial
        other_polynomials += loop_integral.measure.factors

    # insert the regulator to the `other_polynomials` (`numerator` and possily `measure`)
    for poly in other_polynomials:
        poly.polysymbols = poly.polysymbols[:-2] + [loop_integral.regulator] + poly.polysymbols[-2:]
        poly.expolist = np.hstack([poly.expolist[:,:-2], np.zeros([len(poly.expolist),1], dtype=int), poly.expolist[:,-2:]])

    return make_package(
        target_directory = target_directory,
        name = name,

        integration_variables = loop_integral.integration_variables,

        regulators = [loop_integral.regulator],
        requested_orders = [requested_order],

        polynomials_to_decompose = polynomials_to_decompose,
        polynomial_names = ['F','U'],
        other_polynomials = other_polynomials,
        contour_deformation_polynomial = 'F' if contour_deformation else None,

        prefactor = sp.sympify(additional_prefactor) * loop_integral.Gamma_factor * loop_integral.regulator ** loop_integral.regulator_power,

        real_parameters = real_parameters,
        complex_parameters = complex_parameters,

        form_optimization_level = form_optimization_level,
        form_work_space = form_work_space,
        form_insertion_depth = 1, # can have log(UorF(...)) or denominator(UorF(...))

        stabilize = stabilize, 
        decomposition_method = decomposition_method
    )
