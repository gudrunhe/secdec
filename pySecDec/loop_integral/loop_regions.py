import numpy as np
import sympy as sp
import tempfile
from ..misc import sympify_expression
from ..algebra import Polynomial, ExponentiatedPolynomial
from ..polytope import Polytope

def loop_regions(name, loop_integral, smallness_parameter,
                expansion_by_regions_order=0,
                contour_deformation=True,
                additional_prefactor=1, form_optimization_level=2,
                form_work_space='500M',
                form_memory_use=None,
                form_threads=1,
                extra_regulator_name=None,
                extra_regulator_exponent=None,
                decomposition_method='iterative',
                normaliz_executable=None,
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

    :param extra_regulator_name:
        string or sympy symbol, optional;
        Name of the regulator using which monomial factors of the form
        $x_i**(n*p_i)$ are added, to regulate the integrals arising from
        the expansion by regions.

    :param extra_regulator_exponent:
        list of integers or sympy rationals, optional;
        List of the $p_i$ factors of the extra regulator.

    :param decomposition_method:
        string, optional;
        The strategy for decomposing the polynomials. The
        following strategies are available:

        * 'iterative' (default)
        * 'geometric'
        * 'geometric_ku'

    :param normaliz_executable:
        string, optional;
        The command to run `normaliz`. `normaliz` is only
        required if `decomposition_method` starts with
        'geometric'.
        Default: use `normaliz` from pySecDecContrib

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
    polynomials_to_decompose = \
            [loop_integral.exponentiated_U, loop_integral.exponentiated_F] + \
            loop_integral.measure.factors
    
    # add regulators of the form x_i**(n*e_i), where n is the extra regulator
    if extra_regulator_name is not None:
        extra_regulator = sympify_expression(extra_regulator_name)
        loop_integral.regulators.insert(0, extra_regulator)
        if extra_regulator_exponent is not None:
            exponents = extra_regulator_exponent
            if len(exponents) != len(loop_integral.integration_variables):
                raise ValueError(f"extra_regulator_exponent should be a list of {len(loop_integral.integration_variables)} numbers")
        else:
            exponents = suggested_extra_regulator_exponent(
                    loop_integral, smallness_parameter,
                    expansion_by_regions_order=expansion_by_regions_order,
                    normaliz=normaliz_executable)
            print('suggested extra regulator exponents: ', exponents)
        if exponents is not None:
            for exponent, variable in zip(exponents, loop_integral.integration_variables):
                if exponent != 0:
                    variable = Polynomial.from_expression(
                            variable, loop_integral.integration_variables)
                    polynomials_to_decompose.append(ExponentiatedPolynomial(
                        variable.expolist,
                        variable.coeffs,
                        exponent=extra_regulator*exponent,
                        polysymbols=variable.polysymbols
                    ))

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
        numerator = loop_integral.numerator,
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

def suggested_extra_regulator_exponent(loop_integral, smallness_parameter, expansion_by_regions_order=0, normaliz=None):
    '''
    Returns the suggested exponent coefficients for the extra
    regulator needed by the given loop integral, or None if it
    is not needed.

    :param loop_integral:
        :class:`pySecDec.loop_integral`;
        The loop integral which is to be regulated.

    :param smallness_parameter:
        string or sympy symbol;
        The symbol of the variable in which the
        expression is expanded.

    :param expansion_by_regions_order:
        integer;
        The order up to which the expression is expanded
        in the `smallness_parameter`.
        Default: 0

    :param normaliz:
        string, optional;
        The shell command to run `normaliz`.
        Default: use `normaliz` from pySecDecContrib
    '''
    smallness_parameter = sympify_expression(smallness_parameter)
    poly = Polynomial.from_expression(
            str(loop_integral.preliminary_F + loop_integral.preliminary_U),
            loop_integral.preliminary_F.symbols + [smallness_parameter])
    assert poly.symbols[:-1] == loop_integral.integration_variables
    idx = poly.symbols.index(smallness_parameter)
    try:
        constr = extra_regulator_constraints(idx, poly, expansion_by_regions_order,
                loop_integral.regulators, loop_integral.powerlist, loop_integral.dimensionality,
                indices=None, normaliz=normaliz)['all']
        if len(constr) == 0: return None
    except ValueError as e:
        if "need at least one array to concatenate" in str(e):
            return None
        raise e
    assert not constr[:,-1].any()
    constr = constr[:,:-1]
    n = constr.shape[1]
    needed = [constr[:,i].any() for i in range(n)]
    exp = [0] * n
    nzero = n
    while True:
        k = -nzero
        lasti = None
        for i in range(constr.shape[1]):
            if needed[i]:
                exp[i] = 0 if k < 0 else 1 if k < 1 else sp.sympify(1)/sp.prime(k)
                lasti = i
                k += 1
        exp[lasti] = 0
        exp[lasti] = -sum(exp)
        if (constr @ exp).all():
            return exp
        nzero -= 1

def extra_regulator_constraints(exp_param_index, polynomial, exp_order, regulators, powerlist, dimension, indices=None, normaliz=None):
    '''
    Returns a dictionary of vectors :\mathbf{n}_i: that give constraints
    of the form :math:`\langle\mathbf{n_i},\bf{\nu}_{\delta}\rangle \neq 0`
    on the coefficient of a regulator :math:`\delta` introduced by
    multiplying the integrand with a monomial 
    :math:`\mathbf{x}^{\delta\bf{\nu}_{\delta}}`, where :math:`\mathbf{x}`
    are the variables of the input polynomial. The dictionary contains 
    entries for each region individually and a list of all constraints 
    (entry `all`). Only exponents corresponding to integration variables 
    can be non-zero.

    :param exp_param_index:
        int;
        The index of the expansion parameter in the expolist.

    :param polynomial:
        an instance of :class:`.Polynomial`, for which to
        calculate the conditions for.

    :param exp_order:
        int;
        desired expansion order in small parameter; corresponds
        to expansion_by_regions_order in make_regions.

    :param regulators:
        list of regulators (loop_integral.regulators)

    :param powerlist:
        list of propagator exponents (loop_integral.powerlist)

    :param dimension:
        space-time dimension (loop_integral.dimensionality)

    :param indices:
        list of integers or None;
        The indices of the parameters to be included
        in the asymptotic expansion. This should include all
        Feynman parameters (integration variables) and the
        expansion parameter. By default (``indices=None``),
        all parameters are considered.

    :param normaliz:
        string;
        The shell command to run `normaliz`.
        Default: use `normaliz` from pySecDecContrib
    '''

    powerlist = np.array(powerlist)
    powerlist = powerlist[ powerlist > 0 ]
    powerlist = np.insert(powerlist,exp_param_index,0)
    powerlist = np.concatenate((powerlist,[dimension/2]))
    powerlist = sympify_expression(powerlist.tolist())
    for regulator in regulators:
        powerlist = [x.subs(regulator,0) for x in powerlist]

    polytope_vertices = polynomial.expolist

    dim = len(polytope_vertices[0])
    if indices is not None:
        indices = list(indices)
        exp_param_index = indices.index(range(dim)[exp_param_index])
        polytope_vertices = np.array([[vertex[i] for i in indices] for vertex in polytope_vertices])

    polytope = Polytope(vertices=polytope_vertices)
    with tempfile.TemporaryDirectory("normaliz") as tmpdir:
        polytope.complete_representation(normaliz=normaliz, workdir=tmpdir + "/x")

    facets = polytope.facets

    regions = facets[ facets[:,exp_param_index] > 0 ]
    regions = regions[ np.dot(regions,powerlist)<= exp_order ] 

    region_facets_sd_0 = {}
    all_facets_sd_0 = []
    for region in regions:
        # vertices contained in region
        vert = polytope_vertices[ np.inner(polytope_vertices, region[:-1]) + region[-1] == 0 ]
        # project vertices to (small parameter)^0 hyperplane
        vert = np.delete(vert, exp_param_index, 1)
        polytope = Polytope(vertices=vert)
        with tempfile.TemporaryDirectory("normaliz") as tmpdir:
            polytope.complete_representation(normaliz=normaliz, workdir=tmpdir + "/x")
        facets = polytope.facets
        # facets with 0 in aff(facet)
        facets = facets[facets[:,-1]==0][:,:-1]
        region_facets_sd_0[tuple(region[:-1])] = facets
        all_facets_sd_0 = all_facets_sd_0 + [facets]
    all_facets_sd_0 = np.concatenate(all_facets_sd_0)
    
    # internal facets
    facets_sd_int_0 = [facet for facet in all_facets_sd_0  if any([np.array_equal(facet, -fac) for fac in all_facets_sd_0])]
    # store only internal facets for each region
    for region, facets in region_facets_sd_0.items():
        region_facets_sd_0[region] = np.array([facet for facet in facets if any([np.array_equal(facet, fac) for fac in facets_sd_int_0])])
    # select one normal vector for each internal facet
    region_facets_sd_0['all'] = np.array([ facet for facet in facets_sd_int_0 if facet[facet!=0][0]>0 ])

    if indices is not None:
        for region, facets in region_facets_sd_0.items():
            region_facets_sd_0[region] = np.array([[next(facet) if i in indices and i != exp_param_index else 0 for i in range(dim)] for facet in [iter(f) for f in facets]])
    else:
        for region, facets in region_facets_sd_0.items():
            region_facets_sd_0[region] = np.array([[next(facet) if i != exp_param_index else 0 for i in range(dim)] for facet in [iter(f) for f in facets]])

    return region_facets_sd_0
