Overview
========

`pySecDec` is consists of several modules that provide functions and classes for
specific purposes. In this overview, we present only the most important aspects
of selected modules. For detailed instruction of a specific function or class,
please refer to the :doc:`full_reference`.


pySecDec.algebra
----------------

The  `algebra` module implements a very basic computer algebra system. Although
`sympy` in principle provides everything we need, it is way too slow for typical
applications. That is because `sympy` is completely written in `python` without
making use of any precompiled functions. `pySecDec`'s `algebra` module uses the
in general faster `numpy` function wherever possible.

..  _poly_intro:

Polynomials
~~~~~~~~~~~

Since sector decomposition is an algorithm that acts on polynomials, we start with
the key class :class:`Polynomial <pySecDec.algebra.Polynomial>`.
As the name suggests, the :class:`Polynomial <pySecDec.algebra.Polynomial>` class
is a container for multivariate polynomials, i.e. functions of the form:

.. math::

    \sum_i C_i {\prod_j { x_{j}^{\alpha_{ij}} }}

A multivariate polynomial is completely determined by its `coefficients` :math:`C_i`,
the exponents :math:`\alpha_{ij}`. The :class:`Polynomial <pySecDec.algebra.Polynomial>`
class stores these in two arrays:

.. TODO: use doctest

>>> from pySecDec.algebra import Polynomial
>>> poly = Polynomial([[1,0], [0,2]], ['A', 'B'])
>>> poly
 + (A)*x0 + (B)*x1**2
>>> poly.expolist
array([[1, 0],
       [0, 1]])
>>> poly.coeffs
array([A, B], dtype=object)

It is also possible to instatiate the :class:`Polynomial <pySecDec.algebra.Polynomial>`
with by its algebraic representation:

>>> poly2 = Polynomial.from_expression('A*x0 + B*x1**2', ['x0','x1'])
>>> poly2
 + (A)*x0 + (B)*x1**2
>>> poly2.expolist
array([[1, 0],
       [0, 2]])
>>> poly2.coeffs
array([A, B], dtype=object)

Note that the second argument of
:meth:`Polynomial.from_expression() <pySecDec.algebra.Polynomial.from_expression>`
defines the variables :math:`x_j`.

The :class:`Polynomial <pySecDec.algebra.Polynomial>` class implements basic operations:

>>> poly + 1
 + (1) + (B)*x1**2 + (A)*x0
>>> 2 * poly
 + (2*A)*x0 + (2*B)*x1**2
>>> poly + poly
 + (2*B)*x1**2 + (2*A)*x0
>>> poly * poly
 + (B**2)*x1**4 + (2*A*B)*x0*x1**2 + (A**2)*x0**2
>>> poly ** 2
 + (B**2)*x1**4 + (2*A*B)*x0*x1**2 + (A**2)*x0**2


General Expressions
~~~~~~~~~~~~~~~~~~~

In order to perform the :mod:`pySecDec.subtraction` and :pySecDec:`expansion`,
we have to introduce more complex algebraic constructs.

General expressions can be entered in a straightforward way:

>>> from pySecDec.algebra import Expression
>>> log_of_x = Expression('log(x)', ['x'])
>>> log_of_x
log( + (1)*x)

All expressions in the context of this `algebra` module are based
on extending or combining the :class:`Polynomials <pySecDec.algebra.Polynomial>`
introduced :ref:`above <poly_intro>`.
In the example above, ``log_of_x`` is a
:class:`LogOfPolynomial <pySecDec.algebra.LogOfPolynomial>`, which
is a derived class from :class:`Polynomial <pySecDec.algebra.Polynomial>`:

>>> type(log_of_x)
<class 'pySecDec.algebra.LogOfPolynomial'>
>>> isinstance(log_of_x, Polynomial)
True
>>> log_of_x.expolist
array([[1]])
>>> log_of_x.coeffs
array([1], dtype=object)

We have seen an `extension` the
:class:`Polynomial <pySecDec.algebra.Polynomial>` class, now let us consider
a `combination`:

>>> more_complex_expression = log_of_x * log_of_x
>>> more_complex_expression
(log( + (1)*x)) * (log( + (1)*x))

We just introduced the :class:`Product <pySecDec.algebra.Product>`
of two :class:`LogOfPolynomials <pySecDec.algebra.LogOfPolynomial>`:

>>> type(more_complex_expression)
<class 'pySecDec.algebra.Product'>

As suggested before, the :class:`Product <pySecDec.algebra.Product>`
combines two :class:`Polynomials <pySecDec.algebra.Polynomial>`. They
are accessible as the ``factors``:

>>> more_complex_expression.factors[0]
log( + (1)*x)
>>> more_complex_expression.factors[1]
log( + (1)*x)
>>> type(more_complex_expression.factors[0])
<class 'pySecDec.algebra.LogOfPolynomial'>
>>> type(more_complex_expression.factors[1])
<class 'pySecDec.algebra.LogOfPolynomial'>

.. important ::
    When working with this `algebra` module, it is important to understand that
    **everything** is based on the class
    :class:`Polynomials <pySecDec.algebra.Polynomial>`.

To emphasize the importance of the above statement, consider the following code:

>>> expression1 = Expression('x*y', ['x', 'y'])
>>> expression2 = Expression('x*y', ['x'])
>>> type(expression1)
<class 'pySecDec.algebra.Polynomial'>
>>> type(expression2)
<class 'pySecDec.algebra.Polynomial'>
>>> expression1
 + (1)*x*y
>>> expression2
 + (y)*x

Although ``expression1``` and ``expression2`` are mathematically identical,
they are treated differently by the `algebra` module. In ``expression1``, both,
``x`` and ``y``, are considered as variables of the
:class:`Polynomial <pySecDec.algebra.Polynomial>`. In contrast, ``y`` is treated
as `coefficient` in ``expression2``:

>>> expression1.expolist
array([[1, 1]])
>>> expression1.coeffs
array([1], dtype=object)
>>> expression2.expolist
array([[1]])
>>> expression2.coeffs
array([y], dtype=object)

The second argument of the function :func:`Expression <pySecDec.algebra.Expression>`
controls how the variables are distributed between the coefficients and the variables
in the underlying :class:`Polynomials <pySecDec.algebra.Polynomial>`.
Keep that in mind in order to avoid confusion.


Feynman Parametrization of Loop Integrals
-----------------------------------------

The primary purpose of `pySecDec` is calculating loop integrals as they arise in fixed
order calculations in quantum field theories. In our approach, the first step is converting
the loop integral from the momentum representation to the Feynman parameter representation.

.. TODO: give some reference

The module :mod:`pySecDec.loop_integral` implements exactly that conversion.
The most basic use is to calculate the first and second Symanzik polynomials
``U`` and ``F`` from the propagators of a loop integral. Consider for example the
one loop bubble:

.. TODO: check spelling of "Symanzik"

>>> from pySecDec.loop_integral import LoopIntegral
>>> propagators = ['k**2', '(k - p)**2']
>>> loop_momenta = ['k']
>>> one_loop_bubble = LoopIntegral.from_propagators(propagators, loop_momenta)
>>> one_loop_bubble.U
 + (1)*x0 + (1)*x1
>>> one_loop_bubble.F
 + (-p**2)*x0*x1






























