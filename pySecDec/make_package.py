"""
This module provides a make_package like interface to code_writer.sum_package
:func:`.make_sum_package`.

"""

from .code_writer.make_package import MakePackage
from .code_writer.sum_package import sum_package

# ---------------------------------- main function ----------------------------------
def make_package(name, integration_variables, regulators, requested_orders,
                 polynomials_to_decompose, polynomial_names=[], other_polynomials=[],
                 prefactor=1, remainder_expression=1, functions=[], real_parameters=[],
                 complex_parameters=[], form_optimization_level=2, form_work_space='50M',
                 form_memory_use=None, form_threads=2,
                 form_insertion_depth=5, contour_deformation_polynomial=None, positive_polynomials=[],
                 decomposition_method='iterative_no_primary', normaliz_executable='normaliz',
                 enforce_complex=False, split=False, ibp_power_goal=-1, use_iterative_sort=True,
                 use_light_Pak=True, use_dreadnaut=False, use_Pak=True, processes=None, form_executable=None,
                 pylink_qmc_transforms=['korobov3x3']):
    r'''
    Decompose, subtract and expand an expression.
    Return it as c++ package.

    .. seealso::
        In order to decompose a loop integral,
        use the function
        :func:`pySecDec.loop_integral.loop_package`.

    .. seealso::
        The generated library is described in
        :ref:`generated_cpp_libs`.

    .. seealso::
        :func:`pySecDec.code_writer.make_package`

    :param name:
        string;
        The name of the c++ namepace and the output
        directory.

    :param integration_variables:
        iterable of strings or sympy symbols;
        The variables that are to be integrated. The
        intgration region depends on the chosen
        `decomposition_method`.

    :param regulators:
        iterable of strings or sympy symbols;
        The UV/IR regulators of the integral.

    :param requested_orders:
        iterable of integers;
        Compute the expansion in the regulators to these
        orders.

    :param polynomials_to_decompose:
        iterable of strings or sympy expressions or
        :class:`pySecDec.algebra.ExponentiatedPolynomial`
        or :class:`pySecDec.algebra.Polynomial`;
        The polynomials to be decomposed.

    :param polynomial_names:
        iterable of strings;
        Assign symbols for the `polynomials_to_decompose`.
        These can be referenced in the `other_polynomials`;
        see `other_polynomials` for details.

    :param other_polynomials:
        iterable of strings or sympy expressions or
        :class:`pySecDec.algebra.ExponentiatedPolynomial`
        or :class:`pySecDec.algebra.Polynomial`;
        Additional polynomials where no decomposition is
        attempted.
        The symbols defined in `polynomial_names`
        can be used to reference the
        `polynomials_to_decompose`. This is particularly
        useful when computing loop integrals where the
        "numerator" can depend on the first and second
        Symanzik polynomials.

        Example (1-loop bubble with numerator):

        >>> polynomials_to_decompose = ["(x0 + x1)**(2*eps - 4)",
        ...                             "(-p**2*x0*x1)**(-eps))"]
        >>> polynomial_names = ["U", "F"]
        >>> other_polynomials = ["""   (eps - 1)*s*U**2
        ...                          + (eps - 2)*F
        ...                          - (eps - 1)*2*s*x0*U
        ...                          + (eps - 1)*s*x0**2"""]

        .. seealso::
            :mod:`pySecDec.loop_integral`

        Note that the `polynomial_names` refer to the
        `polynomials_to_decompose` **without** their
        exponents.

    :param prefactor:
        string or sympy expression,
        optional;
        A factor that does not depend on the integration
        variables.

    :param remainder_expression:
        string or sympy expression or
        :class:`pySecDec.algebra._Expression`, optional;
        An additional factor.

        Dummy function must be provided with all arguments,
        e.g. ``remainder_expression='exp(eps)*f(x0,x1)'``.
        In addition, all dummy function must be listed
        in `functions`.

    :param functions:
        iterable of strings or sympy symbols, optional;
        Function symbols occuring in `remainder_expression`,
        e.g.``['f']``.

        .. note::
            Only user-defined functions that are provided as
            c++-callable code should be mentioned here.
            Listing basic mathematical functions (e.g. ``log``,
            ``pow``, ``exp``, ``sqrt``, ...) is not required
            and considered an error to avoid name conflicts.

        .. note::
            The power function `pow` and the logarithm `log`
            use the nonstandard continuation with an
            infinitesimal negative imaginary part on the
            negative real axis (e.g. ``log(-1) = -i*pi``).

    :param real_parameters:
        iterable of strings or sympy symbols, optional;
        Symbols to be interpreted as real variables.

    :param complex_parameters:
        iterable of strings or sympy symbols, optional;
        Symbols to be interpreted as complex variables.

    :param form_optimization_level:
        integer out of the interval [0,4], optional;
        The optimization level to be used in FORM.
        Default: ``2``.

    :param form_work_space:
        string, optional;
        The FORM WorkSpace. Default: ``'500M'``.

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

    :param form_insertion_depth:
        nonnegative integer, optional;
        How deep FORM should try to resolve nested function
        calls. Default: ``5``.

    :param contour_deformation_polynomial:
        string or sympy symbol, optional;
        The name of the polynomial in `polynomial_names`
        that is to be continued to the complex plane
        according to a :math:`- i\delta` prescription.
        For loop integrals, this is the second Symanzik
        polynomial ``F``.
        If not provided, no code for contour deformation
        is created.

    :param positive_polynomials:
        iterable of strings or sympy symbols, optional;
        The names of the polynomials in `polynomial_names`
        that should always have a positive real part.
        For loop integrals, this applies to the first Symanzik
        polynomial ``U``.
        If not provided, no polynomial is checked for
        positiveness.
        If `contour_deformation_polynomial` is ``None``, this
        parameter is ignored.

    :param decomposition_method:
        string, optional;
        The strategy to decompose the polynomials. The
        following strategies are available:

        * 'iterative_no_primary' (default): integration region
          :math:`[0,1]^N`.
        * 'geometric_no_primary': integration region :math:`[0,1]^N`.
        * 'geometric_infinity_no_primary': integration region
          :math:`[0,\infty]^N`.
        * 'iterative': primary decomposition followed by
          integration over :math:`[0,1]^{N-1}`.
        * 'geometric': :math:`x_N` is set to one followed by
          integration over :math:`[0,\infty]^{N-1}`.
        * 'geometric_ku': primary decomposition followed by
          integration over :math:`[0,1]^{N-1}`.

        'iterative', 'geometric', and 'geometric_ku' are only
        valid for loop integrals. An end user should use
        'iterative_no_primary', 'geometric_no_primary', or
        'geometric_infinity_no_primary' here.
        In order to compute loop integrals, please use the
        function :func:`pySecDec.loop_integral.loop_package`.

    :param normaliz_executable:
        string, optional;
        The command to run `normaliz`. `normaliz` is only
        required if `decomposition_method` starts with
        'geometric'.
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
        bool or integer, optional;
        Whether or not to split the integration domain
        in order to map singularities from :math:`1` to
        :math:`0`. Set this option to ``True`` if you have
        singularties when one or more integration variables
        are one. If an integer is passed, that integer is
        used as seed to generate the splitting point.
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

    :param form_executable:
        string or None, optional;
        The path to the form exectuable. The argument is passed
        to :meth:`.Coefficient.process`. If ``None``, then either
        ``$FORM``, ``$SECDEC_CONTRIB/bin/form``, or just ``form``
        is used, depending on which environment variable is set.
        Default: ``None``.

    :param pylink_qmc_transforms:
        list or None, optional;
        Required qmc integral transforms, options are:

        * ``korobov<i>x<j>`` for 1 <= i,j <= 6
        * ``korobov<i>`` for 1 <= i <= 6 (same as ``korobov<i>x<i>``)
        * ``sidi<i>`` for 1 <= i <= 6

        `New in version 1.5`.
        Default: ``['korobov3x3']``
    '''

    # Build generators_args
    generators_args = \
    {
        'name' : name + '_integral',
        'integration_variables' : integration_variables,
        'regulators' : regulators,
        'requested_orders' : requested_orders,
        'polynomials_to_decompose' : polynomials_to_decompose,
        'polynomial_names' : polynomial_names,
        'other_polynomials' : other_polynomials,
        'prefactor' : prefactor,
        'remainder_expression' : remainder_expression,
        'functions' : functions,
        'real_parameters' : real_parameters,
        'complex_parameters' : complex_parameters,
        'form_optimization_level' : form_optimization_level,
        'form_work_space' : form_work_space,
        'form_memory_use' : form_memory_use,
        'form_threads' : form_threads,
        'form_insertion_depth' : form_insertion_depth,
        'contour_deformation_polynomial' : contour_deformation_polynomial,
        'positive_polynomials' : positive_polynomials,
        'decomposition_method' : decomposition_method,
        'normaliz_executable' : normaliz_executable,
        'enforce_complex' : enforce_complex,
        'split' : split,
        'ibp_power_goal' : ibp_power_goal,
        'use_iterative_sort' : use_iterative_sort,
        'use_light_Pak' : use_light_Pak,
        'use_dreadnaut' : use_dreadnaut,
        'use_Pak' : use_Pak,
        'processes' : processes,
        'pylink_qmc_transforms' : pylink_qmc_transforms
    }

    sum_package(
                name,
                [MakePackage(**generators_args)],
                regulators,
                requested_orders,
                real_parameters,
                complex_parameters,
                form_executable=form_executable,
                pylink_qmc_transforms=pylink_qmc_transforms
                )
