"""
This module implements a wrapper for loop integrals
around the main program
:func:`pySecDec.code_writer.make_package` - the
function :func:`.loop_package`.

"""

from .from_graph import LoopIntegralFromGraph
from .draw import plot_diagram
from ..code_writer import make_package
from itertools import chain
import numpy as np
import sympy as sp
import os

def loop_package(name, loop_integral, requested_order,
                 real_parameters=[], complex_parameters=[],
                 contour_deformation='auto', force_unsafe_deformation=False,
                 additional_prefactor=1, form_optimization_level=2,
                 form_work_space='500M',
                 decomposition_method='iterative',
                 normaliz_executable='normaliz',
                 normaliz_workdir='normaliz_tmp',
                 enforce_complex=False,
                 split=False):
    '''
    Decompose, subtract and expand a Feynman
    parametrized loop integral. Return it as
    c++ package.

    .. seealso::
        This function is a wrapper around
        :func:`pySecDec.code_writer.make_package`.

    :param name:
        string;
        The name of the c++ namespace and the output
        directory.

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
        bool or 'auto', optional;
        Whether or not to produce code for contour
        deformation. If set to 'auto', produce code
        for contour deformation only if the second
        Symanzik polynomial "F" comes with a negative
        exponent (ignoring the regulator).
        Default: ``'auto'``.

    :param force_unsafe_deformation:
        bool, optional;
        If both Symanzik polynomials "U" and "F" have
        a negative exponent (ignoring the regulator),
        the deformation parameters cannot fully be
        checked for consistency. `pySecDec` always checks
        the :math:`- i\delta` prescription of "F" but it
        cannot check if the contour deformation crosses
        complex poles of "U".
        Default: ``False``.

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

    :param decomposition_method:
        string, optional;
        The strategy to decompose the polynomials. The
        following strategies are available:

        * 'iterative' (default)
        * 'geometric'

        .. note::
            For 'geometric', the third-party program "normaliz"
            is needed. See :ref:`installation_normaliz`.

    :param normaliz_executable:
        string, optional;
        The command to run `normaliz`. `normaliz` is only
        required if `decomposition_method` is set to
        'geometric'.
        Default: 'normaliz'

    :param normaliz_workdir:
        string, optional;
        The working directory for `normaliz`. This directory
        is automatically created and deleted again after
        `normaliz` finishes.
        Default: 'normaliz_tmp'

    :param enforce_complex:
        bool, optional;
        Whether or not the generated integrand functions
        should have a complex return type even though
        they might be purely real.
        The return type of the integrands is automatically
        complex if `contour_deformation` is ``True`` or
        if there are `complex_parameters`. In other cases,
        the calculation can typically be kept purely real.
        Most commonly, this flag is needed if
        ``log(<negative real>)`` occurs in one of the
        integrand functions. However, `pySecDec` will suggest
        setting this flag to ``True`` in that case.
        Default: ``False``

    :param split:
        bool, optional;
        Whether or not to split the integration at :math:`1/2`
        in order to map singularities from :math:`1` to
        :math:`0`. Set this option to ``True`` if you have
        singularties when one or more integration variables
        are one.
        Default: ``False``

    '''
    print('running "loop_package" for "' + name + '"')

    if contour_deformation == 'auto':
        # only need contour deformation of "F" has an overall negative exponent (and a noneuclidean point is calculated)
        contour_deformation = loop_integral.exponent_F.subs(loop_integral.regulator,0) < 0
        print(  'contour deformation is ' + ('required' if contour_deformation else 'not needed')  )
    else:
        contour_deformation = bool(contour_deformation)

    if contour_deformation:
        if loop_integral.exponent_U.subs(loop_integral.regulator,0) < 0: # if "U" in denominator
            if force_unsafe_deformation:
                print(   'WARNING: Contour deformation may silently lead to wrong results by crossing zeros of "U".'   )
            else:
                raise UserWarning('Contour deformation may silently lead to wrong results by crossing zeros of "U". Set ``force_unsafe_deformation=True`` in `loop_package` to proceed anyway.')

    # convert `name` to string
    name = str(name)

    U_and_F = [loop_integral.exponentiated_U.copy(), loop_integral.exponentiated_F.copy()]
    names_U_and_F = sp.symbols(['U','F'])

    # append the regulator symbol and the symbols `U` and `F` to `U` and `F`
    for poly in U_and_F:
        poly.polysymbols.extend([loop_integral.regulator] + names_U_and_F)
        poly.expolist = np.hstack([poly.expolist, np.zeros([len(poly.expolist),len(names_U_and_F)+1], dtype=int)])

    # append the regulator symbol to the `numerator` and to `measure`
    for poly in chain([loop_integral.numerator], loop_integral.measure.factors):
        poly.polysymbols = poly.polysymbols[:-2] + [loop_integral.regulator] + poly.polysymbols[-2:]
        poly.expolist = np.hstack([poly.expolist[:,:-2], np.zeros([len(poly.expolist),1], dtype=int), poly.expolist[:,-2:]])

    other_polynomials = [loop_integral.numerator]

    polynomials_to_decompose = list(U_and_F)
    if sp.sympify( loop_integral.measure ) != 1:
        # need ``loop_integral.measure`` only if it is nontrivial
        polynomials_to_decompose += loop_integral.measure.factors

    make_package_return_value = make_package(
        name = name,

        integration_variables = loop_integral.integration_variables,

        regulators = [loop_integral.regulator],
        requested_orders = [requested_order],

        polynomials_to_decompose = polynomials_to_decompose,
        polynomial_names = names_U_and_F,
        other_polynomials = other_polynomials,
        contour_deformation_polynomial = 'F' if contour_deformation else None,

        prefactor = sp.sympify(additional_prefactor) * loop_integral.Gamma_factor * loop_integral.regulator ** loop_integral.regulator_power,

        real_parameters = real_parameters,
        complex_parameters = complex_parameters,

        form_optimization_level = form_optimization_level,
        form_work_space = form_work_space,

        decomposition_method = decomposition_method,

        normaliz_executable = normaliz_executable,
        normaliz_workdir = normaliz_workdir,

        enforce_complex = enforce_complex,
        split = split
    )

    if isinstance(loop_integral, LoopIntegralFromGraph):
        try:
            plot_diagram(loop_integral.internal_lines, loop_integral.external_lines,
                         os.path.join(name, name), loop_integral.powerlist)
        except Exception as error:
            print('WARNING: Could not draw the Feynman diagram "%s". Reason: %s' % (name,error))

    return make_package_return_value
