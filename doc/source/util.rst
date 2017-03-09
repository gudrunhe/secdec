.. _chapter_secdecutil:

.. tell sphinx that everything is in the namespace `secdecutil`
.. cpp:namespace:: secdecutil

SecDecUtil
==========

`SecDecUtil` is a standalone autotools-c++ package, that collects common helper classes and
functions needed by the c++ code generated using :func:`loop_package <pySecDec.loop_integral.loop_package>`
or :func:`make_package <pySecDec.code_writer.make_package>`. Everything defined by the `SecDecUtil`
is put into the c++ namepace `secdecutil`.



Deep_Apply
----------

Test


Integrators
-----------

Integrand_Container
~~~~~~~~~~~~~~~~~~~

A helper class for containing integrands.

.. cpp:class:: template <typename T, typename ...Args> IntegrandContainer

Integrators/Cuba
~~~~~~~~~~~~~~~~

Integrators/Integrator
~~~~~~~~~~~~~~~~~~~~~~


Pylink
------

Sector_Container
----------------

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
    

Uncertainties
-------------

    .. cpp:class:: template <typename T> UncorrelatedDeviation

        .. cpp:var:: T value

        .. cpp:var:: T uncertainty

This class implements Gaussian uncertainty propagation by overloads of the ``+``, ``-``, ``*``, and ``/`` operators.
It has special overloads for :cpp:class:`std::complex\<T\>`.

.. TODO: Should we add the errors linearly to invalidate the statement below?

Although implemented for completeness, note that division by an :cpp:class:`UncorrelatedDeviation\<std::complex\<T\>\>`
is conceptually wrong, since the denominator is multiplied by its complex conjugate; i.e. two most likely correlated
numbers are multiplied assuming no correlation.
