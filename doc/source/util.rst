.. _chapter_secdecutil:

.. tell sphinx that everything is in the namespace `secdecutil`
.. cpp:namespace:: secdecutil

SecDecUtil
==========

`SecDecUtil` is a standalone autotools-c++ package, that collects common helper classes and
functions needed by the c++ code generated using :func:`loop_package <pySecDec.loop_integral.loop_package>`
or :func:`make_package <pySecDec.code_writer.make_package>`. Everything defined by the `SecDecUtil`
is put into the c++ namespace `secdecutil`.

.. _chapter_secdecutil_amplitude:

Amplitude
---------

A collection of utilities for evaluating amplitudes (sums of integrals multiplied by coefficients).

.. _chapter_secdecutil_weighted_integral:

WeightedIntegral
~~~~~~~~~~~~~~~~

A class template containing an integral, ``I``, and the coefficient of the integral, ``C``. 
A ``WeightedIntegral`` is interpreted as the product ``C*I`` and can be used to represent individual terms in an amplitude.

    .. cpp:struct:: template<typename integral_t, typename coefficient_t> WeightedIntegral

        .. cpp:var:: std::shared_ptr<integral_t> integral;

            A shared pointer to the integral.

        .. cpp:var:: coefficient_t coefficient;

            The coefficient which will be multiplied on to the integral.

        .. cpp:var:: std::string display_name = "WINTEGRAL";

            A string used to indicate the name of the current weighted integral.

        .. cpp:function:: WeightedIntegral(const std::shared_ptr<integral_t>& integral,const coefficient_t& coefficient = coefficient_t(1)) 


The arithmetic operators (``+``, ``-``, ``*``, ``/``) are overloaded for ``WeightedIntegral`` types.

.. _chapter_secdecutil_weighted_integral_handler:

WeightedIntegralHandler
~~~~~~~~~~~~~~~~~~~~~~~

A class template for integrating a sum of ``WeightedIntegral`` types.

    .. cpp:struct:: template<typename integrand_return_t, typename real_t, typename coefficient_t, template<typename...> class container_t> class WeightedIntegralHandler

        .. cpp:var:: bool verbose

            Controls the verbosity of the output of the amplitude.


        .. cpp:var:: real_t min_decrease_factor

            If the next refinement iteration is expected to make the total time taken for the code to run longer than ``wall_clock_limit`` then the number of points to be requested in the next iteration will be reduced by at least ``min_decrease_factor``.

        .. cpp:var:: real_t decrease_to_percentage

            If ``remaining_time * decrease_to_percentage > time_for_next_iteration`` then the number of points requested in the next refinement iteration will be reduced. Here: ``remaining_time = wall_clock_limit - elapsed_time`` and ``time_for_next_iteration`` is the estimated time required for the next refinement iteration. Note: if this condition is met this means that the expected precision will not match the desired precision.

        .. cpp:var:: real_t wall_clock_limit

            If the current elapsed time has passed ``wall_clock`` limit and a refinement iteration finishes then a new refinement iteration will not be started. Instead, the code will return the current result and exit.

        .. cpp:var:: size_t number_of_threads

            The number of threads used to compute integrals concurrently. Note: The integrals themselves may also be computed with multiple threads irrespective of this option.

        .. cpp:var:: size_t reset_cuda_after

            The cuda driver does not automatically remove unnecessary functions from the device memory such that the device may run out of memory after some time. This option controls after how many integrals ``cudaDeviceReset()`` is called to clear the memory. With the default ``0``, ``cudaDeviceReset()`` is never called. This option is ignored if compiled without cuda.

        .. cpp:var:: const container_t<std::vector<term_t>>& expression

            The sum of terms to be integrated.

        .. cpp:var:: real_t epsrel

            The desired relative accuracy for the numerical evaluation of the weighted sum of the sectors.

        .. cpp:var:: real_t epsabs

            The desired absolute accuracy for the numerical evaluation of the weighted sum of the sectors.

        .. cpp:var:: unsigned long long int maxeval

            The maximal number of integrand evaluations for each sector.

        .. cpp:var:: unsigned long long int mineval

            The minimal number of integrand evaluations for each sector.

        .. cpp:var:: real_t maxincreasefac

            The maximum factor by which the number of integrand evaluations will be increased in a single refinement iteration.

        .. cpp:var:: real_t min_epsrel

            The minimum relative accuracy required for each individual sector.

        .. cpp:var:: real_t min_epsabs

            The minimum absolute accuracy required for each individual sector.

        .. cpp:var:: real_t max_epsrel

            The maximum relative accuracy assumed possible for each individual sector. Any sector known to this precision will not be refined further. Note: if this condition is met this means that the expected precision will not match the desired precision.

        .. cpp:var:: real_t max_epsabs

            The maximum absolute accuracy assumed possible for each individual sector. Any sector known to this precision will not be refined further. Note: if this condition is met this means that the expected precision will not match the desired precision.

        .. cpp:var:: ErrorMode errormode
            
            With ``enum ErrorMode : int { abs=0, all, largest, real, imag};``

            Defines how epsrel and epsabs are defined for complex values.
            With the choice  ``largest``, the relative uncertainty is defined as 
            ``max( |Re(error)|, |Im(error)|)/max( |Re(result)|, |Im(result)|)``.
            Choosing ``all`` will apply epsrel and epsabs to both the real
            and imaginary part separately.

.. _chapter_secdecutil_series:

Series
------

A class template for containing (optionally truncated) Laurent series. Multivariate series can be represented as series of series.

This class overloads the arithmetic operators (``+``, ``-``, ``*``, ``/``) and the comparator operators (``==``, ``!=``).
A string representation can be obtained using the ``<<`` operator.
The ``at(i)`` and ``[i]`` operators return the coefficient of the ``i``:sup:`th` power of the expansion parameter. Otherwise elements can be accessed identically to :cpp:class:`std::vector`.

    .. cpp:class:: template <typename T> Series

        .. cpp:var:: std::string expansion_parameter

            A string representing the expansion parameter of the series (default ``x``)

        .. cpp:function:: int get_order_min() const

            Returns the lowest order in the series.

        .. cpp:function:: int get_order_max() const

            Returns the highest order in the series.

        .. cpp:function:: bool get_truncated_above() const

            Checks whether the series is truncated from above.

        .. cpp:function:: bool has_term(int order) const

            Checks whether the series has a term at order ``order``.

        .. cpp:function:: Series(int order_min, int order_max, std::vector<T> content, bool truncated_above = true, const std::string expansion_parameter = "x")


Example:

.. literalinclude:: cpp_doctest/series_doctest_basic_usage.cpp
   :language: cpp

Compile/Run:

.. literalinclude:: cpp_doctest/compile_and_run.txt

Output:

.. literalinclude:: cpp_doctest/series_doctest_basic_usage.txt
   :language: sh

.. _chapter_secdecutil_deep_apply:

Deep Apply
----------

A general concept to apply a :cpp:class:`std::function` to a nested data structure. If the applied :cpp:class:`std::function` is not void then :cpp:func:`deep_apply` returns a nested data structure of the return values. Currently `secdecutil` implements this for :cpp:class:`std::vector` and :cpp:class:`Series`.

This concept allows, for example, the elements of a nested series to be edited without knowing the depth of the nested structure.

    .. cpp:function:: template<typename Tout, typename Tin, template<typename...> class Tnest> Tnest<Tout> deep_apply(Tnest<Tin>& nest, std::function<Tout(Tin)>& func)

Example (complex conjugate a :cpp:class:`Series`):

.. literalinclude:: cpp_doctest/deep_apply_doctest_basic_usage.cpp
   :language: cpp

Compile/Run:

.. literalinclude:: cpp_doctest/compile_and_run.txt

Output:

.. literalinclude:: cpp_doctest/deep_apply_doctest_basic_usage.txt
   :language: sh

.. _chapter_secdecutil_uncertainties:

Uncertainties
-------------

A class template which implements uncertainty propagation for uncorrelated random variables by overloads of the ``+``, ``-``, ``*`` and partially ``/``.
Division by :cpp:class:`UncorrelatedDeviation` is not implemented as it is not always defined. It has special overloads for :cpp:class:`std::complex\<T\>`.

.. note::
    Division by :cpp:class:`UncorrelatedDeviation` is not implemented as this operation is not always well defined.
    Specifically, it is ill defined in the case that the errors are Gaussian distributed as the expectation value,

    .. math::
        \mathrm{E}\left[\frac{1}{X}\right] = \int_{-\infty}^{\infty} \frac{1}{X} p(X)\ \mathrm{d}X,

    where

    .. math::
        p(X) = \frac{1}{\sqrt{2 \pi \sigma^2 }} \exp\left( -\frac{(x-\mu)^2}{2\sigma^2} \right),

    is undefined in the Riemann or Lebesgue sense. The rule :math:`\delta(a/b) = |a/b| \sqrt{ (\delta a/a)^2 + (\delta b/b)^2 }`
    can not be derived from the first principles of probability theory.


The rules implemented for real valued error propagation are:

.. math::
   \delta(a+b) = \sqrt{(\delta a)^2 + (\delta b)^2},

.. math::
   \delta(a-b) = \sqrt{(\delta a)^2 + (\delta b)^2},

.. math::
   \delta(ab) = \sqrt{ (\delta a)^2 b^2 + (\delta b)^2 a^2 + (\delta a)^2 (\delta b)^2 }.

For complex numbers the above rules are implemented for the real and imaginary parts individually.

    .. cpp:class:: template <typename T> UncorrelatedDeviation

        .. cpp:var:: T value

            The expectation value.

        .. cpp:var:: T uncertainty

            The standard deviation.

Example:

.. literalinclude:: cpp_doctest/uncertainties_doctest_basic_usage.cpp
   :language: cpp

Compile/Run:

.. literalinclude:: cpp_doctest/compile_and_run.txt

Output:

.. literalinclude:: cpp_doctest/uncertainties_doctest_basic_usage.txt
   :language: sh

.. _chapter_secdecutil_integrand_container:

Integrand Container
-------------------

A class template for containing integrands. It stores the number of integration variables and the integrand as a :cpp:class:`std::function`.

This class overloads the arithmetic operators (``+``, ``-``, ``*``, ``/``) and the call operator (``()``).

    .. cpp:class:: template <typename T, typename ...Args> IntegrandContainer

        .. cpp:var:: int number_of_integration_variables

            The number of integration variables that the integrand depends on.

        .. cpp:var:: std::function<T(Args...)> integrand

            The integrand function. The call operator forwards to this function.

Example (add two :cpp:class:`IntegrandContainer` and evaluate one point):

.. literalinclude:: cpp_doctest/integrand_container_doctest_basic_usage.cpp
   :language: cpp

Compile/Run:

.. literalinclude:: cpp_doctest/compile_and_run.txt

Output:

.. literalinclude:: cpp_doctest/integrand_container_doctest_basic_usage.txt
   :language: sh

.. _chapter_secdecutil_integrator:

Integrator
----------

A base class template from which integrator implementations inherit. It defines the minimal API available for all integrators.

    .. cpp:class:: template<typename return_t, typename input_t, typename container_t = secdecutil::IntegrandContainer<return_t, input_t const * const>> Integrator

        .. cpp::type:: container_t

            The type of the integrand. It must have the field ``number_of_integration_variables`` and be callable.

        .. cpp:var:: bool together

            (Only available if ``return_t`` is a :cpp:class:`std::complex` type)
            If ``true`` after each call of the function both the real and imaginary parts are passed to the underlying integrator.
            If ``false`` after each call of the function only the real or imaginary part is passed to the underlying integrator.
            For some adaptive integrators considering the real and imaginary part of a complex function separately can improve the sampling.
            Default: ``false``.

        .. cpp:function:: UncorrelatedDeviation<return_t> integrate(const IntegrandContainer<return_t, input_t const * const>&)

            Integrates the :cpp:class:`IntegrandContainer` and returns the value and uncertainty as an :cpp:class:`UncorrelatedDeviation`.

An integrator that chooses another integrator based on the dimension of the integrand.

    .. cpp:class:: template<typename return_t, typename input_t> MultiIntegrator

        .. cpp:var:: Integrator<return_t,input_t>& low_dimensional_integrator

            Reference to the integrator to be used if the integrand has a lower dimension than :cpp:var:`critical_dim`.

        .. cpp:var:: Integrator<return_t,input_t>& high_dimensional_integrator

            Reference to the integrator to be used if the integrand has dimension :cpp:var:`critical_dim` or higher.

        .. cpp:var:: int critical_dim

            The dimension below which the :cpp:var:`low_dimensional_integrator` is used.

.. _chapter_cpp_cquad:

CQuad
~~~~~

For one dimensional integrals, we wrap the cquad integrator form the GNU scientific library (gsl).

CQuad takes the following options:
 * ``epsrel`` -  The desired relative accuracy for the numerical evaluation. Default: ``0.01``.
 * ``epsabs`` - The desired absolute accuracy for the numerical evaluation. Default: ``1e-7``.
 * ``n`` -  The size of the workspace. This value can only be set in the constructor. Changing this attribute of an instance is not possible. Default: ``100``.
 * ``verbose`` -  Whether or not to print status information. Default: ``false``.
 * ``zero_border`` - The minimal value an integration variable can take. Default: ``0.0``. (`new in version 1.3`)

.. _chapter_cpp_qmc:

Qmc
~~~

.. cpp:namespace:: secdecutil::integrators

The quasi-monte carlo integrator as described in [PSD18]_. Using a quasi-monte integrator to compute sector decomposed integrals was pioneered in [LWY+15]_.

.. cpp:class:: template<typename return_t, ::integrators::U maxdim, template<typename,typename,::integrators::U> class transform_t,typename container_t = secdecutil::IntegrandContainer<return_t, typename remove_complex<return_t>::type const * const>,template<typename,typename,::integrators::U> class fitfunction_t = void_template> Qmc : Integrator<return_t,return_t,container_t>, public ::integrators::Qmc<return_t,return_t,maxdim,transform_t,fitfunction_t>

    Derived from :cpp:class:`secdecutil::Integrator` and :cpp:class:`::integrators::Qmc` - the
    underlying standalone implementation of the Qmc.

The most important fields and template arguments of :cpp:class:`Qmc` are:
 * ``minn`` - The minimal number of points in the Qmc lattice. Will be augmented to the next larger available ``n``.
 * ``minm`` - The minimal number of random shifts.
 * ``maxeval`` - The maximal number of integrand evaluations.
 * ``epsrel`` - The desired relative accuracy for the numerical evaluation.
 * ``epsabs`` - The desired absolute accuracy for the numerical evaluation.
 * ``maxdim`` - The highest dimension the :cpp:class:`Qmc` instance can be used for.
 * ``transform_t`` - The periodizing transform to apply prior to integration.
 * ``fitfunction_t`` - The fit function transform to apply for adaptive integration.
 * ``verbosity`` - Controls the amount of status messages during integration. Can be ``0``, ``1``, ``2``, or ``3``.
 * ``devices`` - A :cpp:class:`std::set` of devices to run on. ``-1`` denotes the CPU, positive integers refer to GPUs.

Refer to the documentation of the standalone Qmc for the default values and additional information.

An integral transform has to be chosen by setting the template argument :cpp:type:`transform_t`. Available transforms
are e.g. ``Korobov<r0,r1>`` and ``Sidi<r0>``, please refer to the underlying Qmc implementation for a complete list.
The fit function for adaptive integration can be set by the :cpp:type:`fitfunction_t`, e.g. ``PolySingular``. If not set,
the default of the underlying Qmc implementation is used.

Examples how to use the Qmc :ref:`on the CPU<example_set_qmc_transform_cpp>` and on :ref:`both, CPU and GPU<example_cuda_qmc>` are shown below.

.. cpp:namespace:: secdecutil

.. _chapter_cpp_cuba:

Cuba
~~~~

Currently we wrap the following Cuba integrators:
 * ``Vegas``
 * ``Suave``
 * ``Divonne``
 * ``Cuhre``

The Cuba integrators all implement:
 * ``epsrel`` -  The desired relative accuracy for the numerical evaluation. Default: ``0.01``.
 * ``epsabs`` - The desired absolute accuracy for the numerical evaluation. Default: ``1e-7``.
 * ``flags`` -  Sets the Cuba verbosity flags. The ``flags=2`` means that the Cuba input parameters and the result after each iteration are written to the log file of the numerical integration. Default: ``0``.
 * ``seed`` - The seed used to generate random numbers for the numerical integration with Cuba. Default: ``0``.
 * ``mineval`` -  The number of evaluations which should at least be done before the numerical integrator returns a result. Default: ``0``.
 * ``maxeval`` -  The maximal number of evaluations to be performed by the numerical integrator. Default: ``1000000``.
 * ``zero_border`` - The minimal value an integration variable can take. Default: ``0.0``. (`new in version 1.3`)

The available integrator specific parameters and their default values are:

+----------------------+---------------------+-------------------------+-------------+
| Vegas                | Suave               | Divonne                 | Cuhre       |
+======================+=====================+=========================+=============+
| nstart (``10000``)   | nnew (``1000``)     | key1 (``2000``)         | key (``0``) |
+----------------------+---------------------+-------------------------+-------------+
| nincrease (``5000``) | nmin (``10``)       | key2 (``1``)            |             |
+----------------------+---------------------+-------------------------+-------------+
| nbatch (``500``)     | flatness (``25.0``) | key3 (``1``)            |             |
+----------------------+---------------------+-------------------------+-------------+
|                      |                     | maxpass (``4``)         |             |
+----------------------+---------------------+-------------------------+-------------+
|                      |                     | border (``0.0``)        |             |
+----------------------+---------------------+-------------------------+-------------+
|                      |                     | maxchisq (``1.0``)      |             |
+----------------------+---------------------+-------------------------+-------------+
|                      |                     | mindeviation (``0.15``) |             |
+----------------------+---------------------+-------------------------+-------------+

For the description of these more specific parameters we refer to the Cuba manual.

.. _integrator_examples:

Examples
~~~~~~~~

Integrate Real Function with Cuba Vegas
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Example:

.. literalinclude:: cpp_doctest/integrator_doctest_basic_usage.cpp
   :language: cpp

Compile/Run:

.. literalinclude:: cpp_doctest/compile_and_run_cuba.txt

Output:

.. literalinclude:: cpp_doctest/integrator_doctest_basic_usage.txt
   :language: sh

Integrate Complex Function with Cuba Vegas
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Example:

.. literalinclude:: cpp_doctest/integrator_doctest_complex.cpp
   :language: cpp

Compile/Run:

.. literalinclude:: cpp_doctest/compile_and_run_cuba.txt

Output:

.. literalinclude:: cpp_doctest/integrator_doctest_complex.txt
   :language: sh

Integrate Real Function with Cuba Vegas or CQuad
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Example:

.. literalinclude:: cpp_doctest/integrator_doctest_Vegas_CQuad.cpp
   :language: cpp

Compile/Run:

.. literalinclude:: cpp_doctest/compile_and_run_cuba_gsl.txt

Output:

.. literalinclude:: cpp_doctest/integrator_doctest_Vegas_CQuad.txt
   :language: sh

.. _example_set_qmc_transform_cpp:

Set the integral transform of the Qmc
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Example:

.. literalinclude:: cpp_doctest/integrator_doctest_set_transform_qmc.cpp
   :language: cpp

Compile/Run:

.. literalinclude:: cpp_doctest/compile_and_run_set_transform_qmc.txt

Output:

.. literalinclude:: cpp_doctest/integrator_doctest_set_transform_qmc.txt
   :language: sh

.. _example_cuda_qmc:

Run the Qmc on GPUs
^^^^^^^^^^^^^^^^^^^

Example:

.. literalinclude:: cpp_doctest/integrator_doctest_cuda_qmc.cu
   :language: cpp

Compile/Run:

.. literalinclude:: cpp_doctest/compile_and_run_cuda_qmc.txt

Output:

.. literalinclude:: cpp_doctest/integrator_doctest_cuda_qmc.txt
   :language: sh
