Frequently Asked Questions
==========================

How can I adjust the integrator parameters?
-------------------------------------------

When using the python interface to integration libraries based on *disteval*, the integrator parameters can be specified directly in the integrator call.
Please see the documentation of :class:`DistevalLibrary<pySecDec.integral_interface.DistevalLibrary>` for the list of available parameters.
For example, from ``examples/box1L/integrate_box1L_disteval.py``::

    result = box1L(parameters={"s": 4.0, "t": -0.75, "s1": 1.25, "msq": 1.0},
                   epsrel=1e-3, epsabs=1e-10, format="json")

When using the *disteval* command-line interface, the same parameters can be adjusted via command-line options, as described in :ref:`disteval_cli`.

If the python interface to :class:`IntegralLibrary<pySecDec.integral_interface.IntegralLibrary>` is used for the numerical integration, i.e. a python script like ``examples/box1L/integrate_box1L.py``, the integrator parameters can be specified in the argument list of the integrator call.
For example, using :class:`Qmc<pySecDec.integral_interface.Qmc>` as integrator::

    box1L.use_Qmc(flags=2, epsrel=1e-3, epsabs=1e-12, nstart=5000, nincrease=10000, maxeval=10000000, real_complex_together=True)

The complete list of possible options for the integrators can be found in :mod:`integral_interface <pySecDec.integral_interface>`.

If the C++ interface is used, the options can be specified as fields of the integrator.
For example, after running ``examples/box1L/generate_box1L.py``, in the file ``examples/box1L/integrate_box1L.cpp``, you can modify the corresponding block to e.g.::

    // Integrate
    secdecutil::cuba::Vegas<box1L::integrand_return_t> integrator;
    integrator.flags = 2; // verbose output --> see cuba manual
    integrator.epsrel = 1e-2;
    integrator.epsabs = 1e-12;
    integrator.nstart = 5000;
    integrator.nincrease = 10000;
    integrator.maxeval = 10000000;
    integrator.together = true;

In order to set the Divonne integrator with the same parameters as above, do::

    // Integrate
    secdecutil::cuba::Divonne<box1L::integrand_return_t> integrator;
    integrator.flags = 2; // verbose output --> see cuba manual
    integrator.epsrel = 1e-2;
    integrator.epsabs = 1e-12;
    integrator.maxeval = 10000000;
    integrator.border = 1e-8;
    integrator.together = true;

More information about the C++ integrator class can be found in :numref:`chapter_secdecutil_integrator`.

How can I request a higher numerical accuracy?
----------------------------------------------

The integrator stops if any of the following conditions is fulfilled: (1) ``epsrel`` is reached, (2) ``epsabs`` is reached, (4) timeout is reached (see the ``timeout`` argument in :class:`DistevalLibrary<pySecDec.integral_interface.DistevalLibrary>` and the ``wall_clock_limit`` argument in :class:`IntegralLibrary<pySecDec.integral_interface.IntegralLibrary>`), (4) ``maxeval`` is reached (for :class:`Vegas<pySecDec.integral_interface.Vegas>`, :class:`Suave<pySecDec.integral_interface.Suave>`, :class:`Divonne<pySecDec.integral_interface.Divonne>`, :class:`Cuhre<pySecDec.integral_interface.Cuhre>`, and :class:`Qmc<pySecDec.integral_interface.Qmc>`/:class:`CudaQmc<pySecDec.integral_interface.CudaQmc>`).
Therefore, setting these parameters accordingly will cause the integrator to make more iterations and reach a more accurate result.

What can I do if the integration takes very long?
-------------------------------------------------

For most integrals, the best performance will be achieved using the Quasi-Monte-Carlo integrator from the :ref:`disteval interface<disteval_build>` and we recommend switching to it, if not already used.
If changing the integrator doesn't improve the runtime, it is possible that the integrator parameters should be adjusted, as described in the previous sections.
In particular for integrals with spurious poles, the parameter ``epsabs`` should be increased, since it is the only relevant stopping criterion in this case.

How can I tune the contour deformation parameters?
--------------------------------------------------

Since version 1.5 `pySecDec` automates the selection of contour deformation parameters, and manual tuning is not required in the majority of cases.

Users that wish to experiment can use ``deformation_parameters_maximum``, ``deformation_parameters_minimum``, and ``number_of_presamples`` parameters of :class:`IntegralLibrary <pySecDec.integral_interface.IntegralLibrary>` manually.

What can I do if the program stops with an error message containing `sign_check_error`?
---------------------------------------------------------------------------------------

This error occurs if the contour deformation leads to a wrong sign of the Feynman :math:`i\delta` prescription, usually due to the fact that the deformation parameter :math:`\lambda` is too large.
Since version 1.5 `pySecDec` automates the selection of the :math:`\lambda` parameters, and will automatically adjust the integration contour and retry in cases when `sign_check_error` is detected, so no user intervention is normally required.

If the code still fails to find a valid contour it may display the error message ``All deformation parameters at minimum already, integral still fails`` and stop. In this case try reducing ``deformation_parameters_maximum`` (default: ``1e-5``) to a smaller number. If the code still fails to find a valid contour it may be that your integral has an unavoidable end-point singularity or other numerical problems. Often this error is encountered when the ``real_parameters`` and/or ``complex_parameters`` are very large/small or if some of the parameters differ from each other by orders of magnitude. If all of the ``real_parameters`` or ``complex_parameters`` are of a similar size (but not :math:`\mathcal{O}(1)`) then dividing each parameter by e.g. the largest parameter (such that all parameters are :math:`\mathcal{O}(1)`) can help to avoid a situation where extremely small deformation parameters are required to obtain a valid contour. It may then be possible to restore the desired result using dimensional analysis (i.e. multiplying the result by some power of the largest parameter).

If you still encounter an error after following these suggestions, please open an issue.

What does `additional_prefactor` mean exactly?
----------------------------------------------

We should first point out that the conventions for additional prefactors defined by the user have been changed between `SecDec 3` and `pySecDec`. The prefactor specified by the user will now be *included* in the numerical result.

To make clear what is meant by "additional", we repeat our conventions for Feynman integrals here.

A scalar Feynman graph :math:`G` in :math:`D` dimensions at :math:`L` loops with :math:`N` propagators, where the propagators can have arbitrary, not necessarily integer powers :math:`\nu_j`, has the following representation in momentum space:

.. math::
   :nowrap:

    \begin{align}
    G &= \int\prod\limits_{l=1}^{L} \mathrm{d}^D\kappa_l\;
    \frac{1}
    {\prod\limits_{j=1}^{N} P_{j}^{\nu_j}(\{k\},\{p\},m_j^2)}, \nonumber \\
    \mathrm{d}^D\kappa_l&=\frac{\mu^{4-D}}{i\pi^{\frac{D}{2}}}\,\mathrm{d}^D k_l\;,\;
    P_j(\{k\},\{p\},m_j^2)=(q_j^2-m_j^2+i\delta)\;, \nonumber
    \end{align}

where the :math:`q_j` are linear combinations of external momenta :math:`p_i` and loop momenta :math:`k_l`.

Introducing Feynman parameters leads to:

.. math::

    G = (-1)^{N_{\nu}}
    \frac{\Gamma(N_{\nu}-LD/2)}{\prod_{j=1}^{N}\Gamma(\nu_j)}\int
    \limits_{0}^{\infty}
    \,\prod\limits_{j=1}^{N}dx_j\,\,x_j^{\nu_j-1}\,\delta(1-\sum_{l=1}^N x_l)\,\frac{{\cal U}^{N_{\nu}-(L+1) D/2}}
    {{\cal F}^{N_\nu-L D/2}}

The prefactor :math:`(-1)^{N_{\nu}}\,\Gamma(N_{\nu}-LD/2)/\prod_{j=1}^{N}\Gamma(\nu_j)` coming from the Feynman parametrisation will always be included in the numerical result, corresponding to `additional_prefactor=1` (default), i.e. the program will return the numerical value for :math:`G`. If the user defines `additional_prefactor='gamma(3-2*eps)'`, this prefactor will be expanded in :math:`\epsilon` and included in the numerical result returned by `pySecDec`, in addition to the one coming from the Feynman parametrisation.

For general polynomials not related to loop integrals, i.e. in ``make_package``, the prefactor provided by the user is the only prefactor, as there is no prefactor coming from a Feynman parametrisation in this case. This is the reason why in :func:`make_package <pySecDec.code_writer.make_package>` the keyword for the prefactor defined by the user is ``prefactor``, while in :func:`loop_package <pySecDec.loop_integral.loop_package>` it is ``additional_prefactor``.

.. note::

   The precise normalisation of each output of the python interface is documented :ref:`here <python_interface_prefactor>`.

What can I do if I get `nan`?
-----------------------------

This means that the integral does not converge which can have several reasons. When Divonne is used as an integrator, it is important to use a non-zero value for border, e.g. ``border=1e-8``. Vegas is in general the most robust integrator. When using Vegas, try to increase the values for ``nstart`` and ``nincrease``, for example ``nstart=100000`` (default: ``10000``) and ``nincrease=50000`` (default: ``5000``).

If the integral is non-Euclidean, make sure that `contour_deformation=True` is set.
Another reason for getting `nan` can be that the integral has  singularities at :math:`x_i = 1` and therefore needs usage of the ``split`` option, see item below.

What can I use as numerator of a loop integral?
-----------------------------------------------

The numerator must be a sum of products of numbers, scalar products (e.g. ``p1(mu)*k1(mu)*p1(nu)*k2(nu)`` and/or symbols (e.g. ``m``). The numerator can also be an inverse propagator.
In addition, the numerator must be finite in the limit :math:`\epsilon \rightarrow 0`. The default numerator is ``1``.

Examples::

    p1(mu)*k1(mu)*p1(nu)*k2(nu) + 4*s*eps*k1(mu)*k1(mu)
    p1(mu)*(k1(mu) + k2(mu))*p1(nu)*k2(nu)
    p1(mu)*k1(mu)

More details can be found in :class:`LoopIntegralFromPropagators <pySecDec.loop_integral.LoopIntegralFromPropagators>`.


How can I integrate just one coefficient of a particular order in the regulator?
--------------------------------------------------------------------------------

You can pick a certain order in the C++ interface (see :ref:`cpp_interface`). To integrate only one order, for example the finite part, change the line::

    const box1L::nested_series_t<secdecutil::UncorrelatedDeviation<box1L::integrand_return_t>> result_all = secdecutil::deep_apply( all_sectors, integrator.integrate );

to::

    int order = 0; // compute finite part only
    const secdecutil::UncorrelatedDeviation<box1L::integrand_return_t> result_order = secdecutil::deep_apply(all_sectors.at(order), integrator.integrate );

where ``box1L`` is to be replaced by the name of your integral. In addition, you should change the lines::

    std::cout << "-- integral without prefactor -- " << std::endl;
    std::cout << result_all << std::endl << std::endl;

to::

    std::cout << "-- integral without prefactor -- " << std::endl;
    std::cout << result_order << std::endl << std::endl;

and remove the lines::

    std::cout << "-- prefactor -- " << std::endl;
    const box1L::nested_series_t<box1L::integrand_return_t> prefactor = box1L::prefactor(real_parameters, complex_parameters);
    std::cout << prefactor << std::endl << std::endl;

    std::cout << "-- full result (prefactor*integral) -- " << std::endl;
    std::cout << prefactor*result_all << std::endl;

because the expansion of the prefactor will in general mix with the pole coefficients and thus affect the finite part. We should point out however that deleting these lines also means that the result will not contain any prefactor, not even the one coming from the Feynman parametrisation.

How can I use complex masses?
-----------------------------

In the python script generating the expressions for the integral, define mass symbols in the same way as for real masses, e.g::

    Mandelstam_symbols=['s']
    mass_symbols=['msq']

Then, in :mod:`loop_package <pySecDec.loop_integral.loop_package>` define::

    real_parameters = Mandelstam_symbols,
    complex_parameters = mass_symbols,

In the integration script (using the python interface), the numerical values for the complex parameters are given after the ones for the real parameters::

    str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = integral(real_parameters=[4.],complex_parameters=[1.-0.0038j])

Note that in python the letter ``j`` is used rather than ``i`` for the imaginary part.

In the C++ interface, you can set (for the example `triangle2L`)::

    const std::vector<triangle2L::real_t> real_parameters = { 4. };
    const std::vector<triangle2L::complex_t> complex_parameters = { {1.,0.0038} };


When should I use the "split" option?
-------------------------------------

The modules :func:`loop_package <pySecDec.loop_integral.loop_package>` and :func:`make_package <pySecDec.code_writer.make_package>` have the option to split the integration domain (``split=True``). This option can be useful for integrals which do not have a Euclidean region. If certain kinematic conditions are fulfilled, for example if the integral contains massive on-shell lines, it can happen that singularities at :math:`x_i = 1` remain in the :math:`\mathcal{F}` polynomial after the decomposition. The split option remaps these singularities to the origin of parameter space. If your integral is of this type, and with the standard approach the numerical integration does not seem to converge, try the ``split`` option. It produces a lot more sectors, so it should not be used without need. We also would like to mention that very often a change of basis to increase the (negative) power of the :math:`\mathcal{F}` polynomial can be beneficial if integrals of this type occur in the calculation.

How can I obtain results from pySecDec in a format convenient for GiNaC/ Sympy/ Mathematica/ Maple?
---------------------------------------------------------------------------------------------------

When using the python interface to *disteval* libraries (i.e. :class:`DistevalLibrary<pySecDec.integral_interface.DistevalLibrary>`), use the ``format`` argument to change the output format. Similar format selection is available in the *disteval* command-line interface, as described in :ref:`disteval_cli`.

When using :class:`IntegralLibrary<pySecDec.integral_interface.IntegralLibrary>`), you can use the functions :func:`series_to_ginac <pySecDec.integral_interface.series_to_ginac>`, :func:`series_to_sympy <pySecDec.integral_interface.series_to_sympy>`, :func:`series_to_mathematica <pySecDec.integral_interface.series_to_mathematica>`, :func:`series_to_maple <pySecDec.integral_interface.series_to_maple>` to convert the output as needed.

Example::

    #!/usr/bin/env python3
    from pySecDec.integral_interface import IntegralLibrary
    from pySecDec.integral_interface import series_to_ginac, series_to_sympy, series_to_mathematica, series_to_maple

    if __name__ == "__main__":

        # load c++ library
        easy = IntegralLibrary('easy/easy_pylink.so')

        # integrate
        _, _, result = easy()

        # print result
        print(series_to_ginac(result))
        print(series_to_sympy(result))
        print(series_to_mathematica(result))
        print(series_to_maple(result))

Outputs::

    ('(1+0*I)/eps + (0.306852819440052549+0*I) + Order(eps)', '(5.41537065611170534e-17+0*I)/eps + (1.3864926114078559e-15+0*I) + Order(eps)')
    ('(1+0*I)/eps + (0.306852819440052549+0*I) + O(eps)', '(5.41537065611170534e-17+0*I)/eps + (1.3864926114078559e-15+0*I) + O(eps)')
    ('(1+0*I)/eps + (0.306852819440052549+0*I) + O[eps]', '(5.41537065611170534*10^-17+0*I)/eps + (1.3864926114078559*10^-15+0*I) + O[eps]')
    ('(1+0*I)/eps + (0.306852819440052549+0*I) + O(eps)', '(5.41537065611170534e-17+0*I)/eps + (1.3864926114078559e-15+0*I) + O(eps)')


Expansion by regions: what does the parameter ``z`` mean?
---------------------------------------------------------

When expansion by regions via the "rescaling with z-method" is used, the parameter ``z`` acts as expansion parameter in the Taylor expansion of the integrand. After the code generation step, in the numerical integration, ``z=1`` needs to be used and the kinematic invariants have to be set to the same values as would be used with the t-method, i.e. the kinematic values desired by the user.

Expansion by regions: why does the t-method not converge?
---------------------------------------------------------

With the t-method, configurations can occur for particular kinematic points which, after sector decomposition, lead to a pole at the upper integration boundary, where the contour deformation vanishes and therefore cannot regulate this pole.
In such a case the z-method should be used, because it does not transform the Feynman parameters in a way which can induce such a configuration.

.. _python_interface_prefactor:

What exactly is returned when calling :func:`IntegralLibrary <pySecDec.integral_interface.IntegralLibrary>` in the python interface?
------------------------------------------------------------------------------------------------------------------------------------

In order to compute an integral in the python interface, first an :func:`IntegralLibrary <pySecDec.integral_interface.IntegralLibrary>` needs to be instantiated and then called.

.. code::

    myintegral = IntegralLibrary('myintegral/myintegral_pylink.so')
    str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = myintegral()

The call to ``myintegral`` will return 3 strings. 
The precise definition each string depends on how the integral library was initially generated.
In the above code, ``str_prefactor`` returns ``1`` in *almost* all cases.
The one exception is if the integral library was generated by a call to the ``code_writer`` version of ``make_package``.
If you are using any of the primarily user facing functions (i.e. those we demonstrate in the ``examples/`` folder), then we would advise just using ``str_integral_with_prefactor`` and not utilising ``str_integral_without_prefactor`` or ``str_prefactor``.

Below we document the various ways an integral library can be generated and the precise content of the prefactor string.
Note that any internally generated prefactors, for example from Feynman parametrisation, are *always* included in the integral and not ``str_prefactor``.

:func:`pySecDec.sum_package <pySecDec.code_writer.sum_package>` (:func:`pySecDec.code_writer.sum_package`)

    The prefactor may be specified in the ``package_generator`` passed to :func:`pySecDec.sum_package <pySecDec.code_writer.sum_package>`, the two options are:

    * ``MakePackage(prefactor='x', ...)`` :  ``x`` is included in ``str_integral_without_prefactor`` and ``str_prefactor = 1``

    * ``LoopPackage(additional_prefactor='x', ...)`` : ``x`` is included in ``str_integral_without_prefactor`` and ``str_prefactor = 1``

:func:`pySecDec.make_package` 

    This function is a thin wrapper around :func:`pySecDec.code_writer.sum_package`

    * ``make_package(prefactor='x', ...)`` : ``x`` is included in ``str_integral_without_prefactor`` and ``str_prefactor = 1``

:func:`pySecDec.loop_package <pySecDec.loop_integral.loop_package>` (:func:`pySecDec.loop_integral.loop_package`)

    This function is a thin wrapper around :func:`pySecDec.make_package`

    * ``loop_package(additional_prefactor='x', ...)`` : ``x`` is included in ``str_integral_without_prefactor`` and ``str_prefactor = 1``

:func:`pySecDec.code_writer.make_package`

    This function is primarily used internally to handle the generation of sector decomposed integral libraries.
    
    * ``make_package(prefactor='x', ...)`` : ``x`` is not included in ``str_integral_without_prefactor`` and ``str_prefactor = x``, ``str_integral_with_prefactor = x * str_integral_without_prefactor``.

What should I do if I get ``OSError: [Errno 24] Too many open files``?
---------------------------------------------------------

This happens if the number of running subprocesses exceeds the open file limit of the OS.
Due to the massive parallelization during integration, this can happen if the default limit is set too low.
How to raise the open file limit depends on the OS. On Mac OSX (El Capitan), where this is known to be a problem, 
the command #ulimit -Sn 10000 sets the limit to 10000.

