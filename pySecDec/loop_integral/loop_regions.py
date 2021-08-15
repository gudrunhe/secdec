import numpy as np
import sympy as sp
from ..misc import sympify_expression
from pySecDec.algebra import Polynomial, ExponentiatedPolynomial

def loop_regions(name, loop_integral, smallness_parameter,
                expansion_by_regions_order=0,
                contour_deformation=True,
                additional_prefactor=1, form_optimization_level=2,
                form_work_space='500M',
                add_monomial_regulator_power=None,
                decomposition_method='iterative',
                normaliz_executable='normaliz',
                enforce_complex=False,
                split=False, ibp_power_goal=-1,
                use_iterative_sort=True, use_light_Pak=True,
                use_dreadnaut=False, use_Pak=True,
                processes=None):
    """
    Apply expansion by regions method to the loop integral using the function.

    .. seealso::
        This function is a wrapper around
        :func:`pySecDec.make_regions`.

    .. seealso::
        The generated library is described in
        :ref:`generated_cpp_libs`.

    :param name:
        string;
        The name of the c++ namespace and the output
        directory.

    :param loop_integral:
        :class:`pySecDec.loop_integral`;
        The loop integral to which the expansion by regions method is applied.

    :param smallness_parameter:
        string or sympy symbol;
        The symbol of the variable in which the
        expression is expanded.

    :param expansion_by_regions_order:
        integer;
        The order up to which the expression is expanded
        in the `smallness_parameter`.
        Default: 0

    :param contour_deformation:
        bool, optional;
        Whether or not to produce code for contour
        deformation.
        Default: ``True``.

    :param additional_prefactor:
        string or sympy expression, optional;
        An additional factor to be multiplied to the loop
        integral. It may depend on the regulators, the
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

    :param add_monomial_regulator_power:
        string or sympy symbol;
        Name of the regulator, using which monomial factors of the form
        $x_i**(n/p_i)$ are added, to regulate the integrals arising from
        the expansion by regions.

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


    """
    polynomials_to_decompose = [loop_integral.exponentiated_U, loop_integral.exponentiated_F] + loop_integral.measure.factors + [loop_integral.numerator]
    
    # add regulators of the form x_i**(n/p_i), where n is a regulator
    if add_monomial_regulator_power is not None:
        regulator = sympify_expression(add_monomial_regulator_power)
        loop_integral.regulators.insert(0,regulator)
        primes = [sp.prime(n+1) for n in range(len(loop_integral.integration_variables)-1)]
        monomial_factors = []
        for prime, variable in zip(primes,loop_integral.integration_variables[:-1]):
            variable = Polynomial.from_expression(variable,loop_integral.integration_variables)
            monomial = ExponentiatedPolynomial(variable.expolist,variable.coeffs,exponent=regulator/prime, polysymbols=variable.polysymbols)
            monomial_factors.append(monomial)
        variable = Polynomial.from_expression(loop_integral.integration_variables[-1],loop_integral.integration_variables)
        monomial_factors.append(ExponentiatedPolynomial(variable.expolist,variable.coeffs,exponent=sum(-regulator/prime for prime in primes), polysymbols=variable.polysymbols))
        polynomials_to_decompose = [loop_integral.exponentiated_U, loop_integral.exponentiated_F] + loop_integral.measure.factors + monomial_factors + [loop_integral.numerator]

    polynomial_names = ["U","F"] + ['Poly%i' % i for i in range(len(polynomials_to_decompose)-2)]

    if "make_regions" not in dir():
        from ..make_regions import make_regions

    return make_regions(
        # make_regions_args
        name = name,
        integration_variables = loop_integral.integration_variables,
        regulators = loop_integral.regulators,
        requested_orders = [],
        smallness_parameter = smallness_parameter,
        polynomials_to_decompose = polynomials_to_decompose,
        expansion_by_regions_order = expansion_by_regions_order,
        real_parameters = [],
        complex_parameters = [],
        polytope_from_sum_of = [0,1],
        normaliz=normaliz_executable,

        # make_package_args
        polynomial_names = polynomial_names,
        prefactor = sympify_expression(additional_prefactor) * loop_integral.Gamma_factor,

        contour_deformation_polynomial = 'F' if contour_deformation else None,
        positive_polynomials = ['U'],

        form_optimization_level = form_optimization_level,
        form_work_space = form_work_space,

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