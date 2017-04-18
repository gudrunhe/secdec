.. _chapter_secdecutil:

.. tell sphinx that everything is in the namespace `secdecutil`
.. cpp:namespace:: secdecutil

SecDecUtil
==========

`SecDecUtil` is a standalone autotools-c++ package, that collects common helper classes and
functions needed by the c++ code generated using :func:`loop_package <pySecDec.loop_integral.loop_package>`
or :func:`make_package <pySecDec.code_writer.make_package>`. Everything defined by the `SecDecUtil`
is put into the c++ namepace `secdecutil`.

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

Integrand Container
-------------------

A class template for containing integrands. It stores the number of integration variables and the integrand as a :cpp:class:`std::function`.

This class overloads the arithmetic operators (``+``, ``-``, ``*``, ``/``).

    .. cpp:class:: template <typename T, typename ...Args> IntegrandContainer

        .. cpp:var:: int number_of_integration_variables
      
            The number of integration variables that the integrand depends on.

        .. cpp:var:: std::function<T(Args...)> integrand

            The integrand function.

Example (add two :cpp:class:`IntegrandContainer` and evaluate one point):

.. literalinclude:: cpp_doctest/integrand_container_doctest_basic_usage.cpp
   :language: cpp

Compile/Run: 

.. literalinclude:: cpp_doctest/compile_and_run.txt
   
Output:

.. literalinclude:: cpp_doctest/integrand_container_doctest_basic_usage.txt
   :language: sh

Integrator
----------

A base class template from which integrator implementations inherit. It defines the minimal API available for all integrators.

    .. cpp:class:: template<typename return_t, typename input_t> Integrator

        .. cpp:var:: bool together 
      
            (Only available if ``return_t`` is a :cpp:class:`std::complex` type) 
            If ``true`` after each call of the function both the real and imaginary parts are passed to the underlying integrator.
            If ``false`` after each call of the function only the real or imaginary part is passed to the underlying integrator.
            For some adaptive integrators considering the real and imaginary part of a complex function separately can improve the sampling.
            Default: ``false``.

        .. cpp:function:: UncorrelatedDeviation<return_t> integrate(const IntegrandContainer<return_t, input_t const * const>&)

            Integrates the :cpp:class:`IntegrandContainer` and returns the value and uncertainty as an :cpp:class:`UncorrelatedDeviation`.

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

The available integrator specific parameters and their default values are:

+---------------------+---------------------+-------------------------+-------------+
| Vegas               | Suave               | Divonne                 | Cuhre       |
+=====================+=====================+=========================+=============+
| nstart (``1000``)   | nnew (``1000``)     | key1 (``2000``)         | key (``0``) |
+---------------------+---------------------+-------------------------+-------------+
| nincrease (``500``) | nmin (``10``)       | key2 (``1``)            |             |
+---------------------+---------------------+-------------------------+-------------+
| nbatch (``500``)    | flatness (``25.0``) | key3 (``1``)            |             |
+---------------------+---------------------+-------------------------+-------------+
|                     |                     | maxpass (``4``)         |             |
+---------------------+---------------------+-------------------------+-------------+
|                     |                     | border (``0.0``)        |             |
+---------------------+---------------------+-------------------------+-------------+
|                     |                     | maxchisq (``1.0``)      |             |
+---------------------+---------------------+-------------------------+-------------+
|                     |                     | mindeviation (``0.15``) |             |
+---------------------+---------------------+-------------------------+-------------+

For the description of these more specific parameters we refer to the Cuba manual.

Examples
~~~~~~~~

Integrate Real Function with Cuba Vegas
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Example:

.. literalinclude:: cpp_doctest/integrator_doctest_basic_usage.cpp
   :language: cpp

Compile/Run: 

.. literalinclude:: cpp_doctest/compile_and_run_Cuba.txt
   
Output:

.. literalinclude:: cpp_doctest/integrator_doctest_basic_usage.txt
   :language: sh

Integrate Complex Function with Cuba Vegas
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Example:

.. literalinclude:: cpp_doctest/integrator_doctest_complex.cpp
   :language: cpp

Compile/Run: 

.. literalinclude:: cpp_doctest/compile_and_run_Cuba.txt
   
Output:

.. literalinclude:: cpp_doctest/integrator_doctest_complex.txt
   :language: sh
