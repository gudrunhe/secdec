"""
This module implements a wrapper for loop integrals
around the main program
:func:`pySecDec.code_writer.make_package` - the
function :func:`.loop_package`.

"""

from .from_graph import LoopIntegralFromGraph
from .draw import plot_diagram
from ..algebra import Polynomial
from ..code_writer import make_package
from ..misc import flatten, sympify_expression
from itertools import chain
import numpy as np
import sympy as sp
import os

def loop_package(name, loop_integral, requested_order,
                 real_parameters=[], complex_parameters=[],
                 contour_deformation=True,
                 additional_prefactor=1, form_optimization_level=2,
                 form_work_space='50M',
                 form_memory_use=None,
                 form_threads=2,
                 decomposition_method='iterative',
                 normaliz_executable='normaliz',
                 enforce_complex=False,
                 split=False, ibp_power_goal=-1,
                 use_iterative_sort=True, use_light_Pak=True,
                 use_dreadnaut=False, use_Pak=True,
                 processes=None):
    '''
    Decompose, subtract and expand a Feynman
    parametrized loop integral. Return it as
    c++ package.

    .. seealso::
        This function is a wrapper around
        :func:`pySecDec.code_writer.make_package`.

    .. seealso::
        The generated library is described in
        :ref:`generated_cpp_libs`.

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

    :param complex_parameters:
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
        string or sympy expression, optional;
        An additional factor to be multiplied to the loop
        integral. It may depend on the regulator, the
        `real_parameters`, and the `complex_parameters`.

    :param form_optimization_level:
        integer out of the interval [0,4], optional;
        The optimization level to be used in FORM.
        Default: ``2``.

    :param form_work_space:
        string, optional;
        The FORM WorkSpace. Default: ``'50M'``.

        Setting this to smaller values will reduce FORM memory
        usage (without affecting performance), but each problem
        has some minimum value below which FORM will refuse to
        work: it will fail with error message indicating that
        larger WorkSpace is needed, at which point WorkSpace
        will be adjusted and FORM will be re-run.

    :param form_memory_use:
        string, optional;
        The target FORM memory usage. When specified, `form.set`
        parameters will be adjusted so that FORM uses at most
        approximately this much resident memory.

        The minimum is approximately to 600M + 350M per worker
        thread if ``form_work_space`` is left at ``'50M'``.
        if form_work_space is increased to ``'500M'``, then
        the minimum is 2.5G + 2.5G per worker thread.
        Default: ``None``, meaning use the default FORM values.

    :param form_threads:
        integer, optional;
        Number of threads (T)FORM will use. Default: ``2``.

    :param decomposition_method:
        string, optional;
        The strategy for decomposing the polynomials. The
        following strategies are available:

        * 'iterative' (default)
        * 'geometric'
        * 'geometric_ku'

        .. note::
            For 'geometric' and 'geometric_ku', the
            third-party program "normaliz" is needed.
            See :ref:`installation_normaliz`.

    :param normaliz_executable:
        string, optional;
        The command to run `normaliz`. `normaliz` is only
        required if `decomposition_method` is set to
        'geometric' or 'geometric_ku'.
        Default: 'normaliz'

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
        Whether or not to split the integration domain in
        order to map singularities from :math:`1` to
        :math:`0`. Set this option to ``True`` if you have
        singularties when one or more integration variables
        are one.
        Default: ``False``

    :param ibp_power_goal:
        number or iterable of number, optional;
        The `power_goal` that is forwarded to
        :func:`.integrate_by_parts`.

        This option controls how the subtraction terms are
        generated. Setting it to ``-numpy.inf`` disables
        :func:`.integrate_by_parts`, while ``0`` disables
        :func:`.integrate_pole_part`.

        .. versionadded: 1.4
            A separate power_goal for each of the
            `integration_variables` can be set by passing an
            iterable.

        .. seealso::
            To generate the subtraction terms, this function
            first calls :func:`.integrate_by_parts` for each
            integration variable with the give `ibp_power_goal`.
            Then :func:`.integrate_pole_part` is called.

        Default: ``-1``

    :param use_iterative_sort:
        bool;
        Whether or not to use
        :func:`.squash_symmetry_redundant_sectors_sort`
        with :func:`.iterative_sort` to find sector symmetries.
        Default: ``True``

    :param use_light_Pak:
        bool;
        Whether or not to use
        :func:`.squash_symmetry_redundant_sectors_sort`
        with :func:`.light_Pak_sort` to find sector symmetries.
        Default: ``True``

    :param use_dreadnaut:
        bool or string, optional;
        Whether or not to use
        :func:`.squash_symmetry_redundant_sectors_dreadnaut`
        to find sector symmetries.
        If given a string, interpret that string as the command
        line executable `dreadnaut`. If ``True``, try
        ``$SECDEC_CONTRIB/bin/dreadnaut`` and, if the
        environment variable ``$SECDEC_CONTRIB`` is not set,
        ``dreadnaut``.
        Default: ``False``

    :param use_Pak:
        bool;
        Whether or not to use
        :func:`.squash_symmetry_redundant_sectors_sort`
        with :func:`.Pak_sort` to find sector symmetries.
        Default: ``True``

    :param processes:
        integer or None, optional;
        The maximal number of processes to be used. If ``None``,
        the number of CPUs :func:`multiprocessing.cpu_count()` is
        used.
        `New in version 1.3`.
        Default: ``None``

    '''
    print('running "loop_package" for "' + name + '"')

    # convert `contour_deformation` to bool
    contour_deformation = bool(contour_deformation)

    # convert `name` to string
    name = str(name)

    U_and_F = [loop_integral.exponentiated_U.copy(), loop_integral.exponentiated_F.copy()]
    names_U_and_F = sp.symbols(['U','F'])

    # append the regulator symbol and the symbols `U` and `F` to `U` and `F`
    for poly in U_and_F:
        poly.polysymbols.extend([loop_integral.regulator] + names_U_and_F)
        poly.expolist = np.hstack([poly.expolist, np.zeros([len(poly.expolist),len(names_U_and_F)+1], dtype=int)])

    # append the regulator symbol to the `numerator` and to `measure`
    numerator = loop_integral.numerator.copy()
    measure = loop_integral.measure.copy()
    for poly in chain([numerator], measure.factors):
        poly.polysymbols = poly.polysymbols[:-2] + [loop_integral.regulator] + poly.polysymbols[-2:]
        poly.expolist = np.hstack([poly.expolist[:,:-2], np.zeros([len(poly.expolist),1], dtype=int), poly.expolist[:,-2:]])

    if np.issubdtype(numerator.coeffs.dtype, np.number):
        other_polynomials = [numerator]
    else:
        symbols = numerator.polysymbols
        numerator.coeffs = np.array( [Polynomial.from_expression(coeff, symbols) for coeff in numerator.coeffs] )
        other_polynomials = [flatten(numerator, 1)]

    polynomials_to_decompose = list(U_and_F)
    if sympify_expression( measure ) != 1:
        # need ``loop_integral.measure`` only if it is nontrivial
        polynomials_to_decompose += measure.factors

    make_package_return_value = make_package(
        name = name,

        integration_variables = loop_integral.integration_variables,

        regulators = [loop_integral.regulator],
        requested_orders = [requested_order],

        polynomials_to_decompose = polynomials_to_decompose,
        polynomial_names = names_U_and_F,
        other_polynomials = other_polynomials,
        contour_deformation_polynomial = 'F' if contour_deformation else None,
        positive_polynomials = ['U'],

        prefactor = sympify_expression(additional_prefactor) * loop_integral.Gamma_factor * loop_integral.regulator ** loop_integral.regulator_power,

        real_parameters = real_parameters,
        complex_parameters = complex_parameters,

        form_optimization_level = form_optimization_level,
        form_work_space = form_work_space,
        form_memory_use = form_memory_use,
        form_threads = form_threads,

        decomposition_method = decomposition_method,

        normaliz_executable = normaliz_executable,

        use_iterative_sort = use_iterative_sort,
        use_Pak = use_Pak,
        use_dreadnaut = use_dreadnaut,
        use_light_Pak = use_light_Pak,

        enforce_complex = enforce_complex,
        ibp_power_goal = ibp_power_goal,
        split = split,
        processes = processes
    )

    if isinstance(loop_integral, LoopIntegralFromGraph):
        try:
            plot_diagram(loop_integral.internal_lines, loop_integral.external_lines,
                         os.path.join(name, name), loop_integral.powerlist)
        except Exception as error:
            print('WARNING: Could not draw the Feynman diagram "%s". Reason: %s' % (name,error))

    return make_package_return_value
