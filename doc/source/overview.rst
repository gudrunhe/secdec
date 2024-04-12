Overview
========

`pySecDec` consists of several modules that provide functions and classes for
specific purposes. In this overview, we present only the most important aspects
of selected modules. These are exactly the modules necessary to set up the algebraic 
computation of a Feynman loop integral requisite for the numerical evaluation. 
For detailed instruction of a specific function or class,
please be referred to the :doc:`reference guide <full_reference>`.


.. _algebra_intro:

The Algebra Module
------------------

The  `algebra` module implements a very basic computer algebra system. `pySecDec` 
uses both `sympy` and `numpy`. Although
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

A multivariate polynomial is completely determined by its `coefficients` :math:`C_i` and 
the exponents :math:`\alpha_{ij}`. The :class:`Polynomial <pySecDec.algebra.Polynomial>`
class stores these in two arrays:

>>> from pySecDec.algebra import Polynomial
>>> poly = Polynomial([[1,0], [0,2]], ['A', 'B'])
>>> poly
 + (A)*x0 + (B)*x1**2
>>> poly.expolist
array([[1, 0],
       [0, 2]])
>>> poly.coeffs
array([A, B], dtype=object)

It is also possible to instantiate the :class:`Polynomial <pySecDec.algebra.Polynomial>`
by its algebraic representation:

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

Within the :class:`Polynomial <pySecDec.algebra.Polynomial>` class, basic operations are implemented:

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

In order to perform the :mod:`pySecDec.subtraction` and :mod:`pySecDec.expansion`,
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

We have seen an `extension` to the
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
are accessible through the ``factors``:

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
    :class:`Polynomial <pySecDec.algebra.Polynomial>`.

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

Although ``expression1`` and ``expression2`` are mathematically identical,
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
controls how the variables are distributed among the coefficients and the variables
in the underlying :class:`Polynomials <pySecDec.algebra.Polynomial>`.
Keep that in mind in order to avoid confusion. One can always check which symbols are
considered as variables by asking for the ``symbols``:

>>> expression1.symbols
[x, y]
>>> expression2.symbols
[x]


.. _loop_integral:

Feynman Parametrization of Loop Integrals
-----------------------------------------

The primary purpose of `pySecDec` is the numerical computation of loop integrals as they arise in fixed
order calculations in quantum field theories. 

The conventions of `pySecDec` are fixed as follows:
a Feynman graph :math:`G^{\mu_1 \ldots \mu_R}_{l_1 \ldots l_R}` in :math:`D` dimensions at :math:`L` loops with :math:`R` loop momenta in the numerator and :math:`N` propagators, where the propagators can have arbitrary, not necessarily integer powers :math:`\nu_j`, is considered to have the following representation in momentum space,

.. math::
   :nowrap:

    \begin{align}
    G &= \int\prod\limits_{l=1}^{L} \mathrm{d}^D\kappa_l\;
    \frac{k_{l_1}^{\mu_1} \cdots k_{l_R}^{\mu_R}}
    {\prod\limits_{j=1}^{N} P_{j}^{\nu_j}(\{k\},\{p\},m_j^2)}, \nonumber \\
    \mathrm{d}^D\kappa_l&=\frac{\mu^{4-D}}{i\pi^{\frac{D}{2}}}\,\mathrm{d}^D k_l\;,\;
    P_j(\{k\},\{p\},m_j^2)=(q_j^2-m_j^2+i\delta)\;, \nonumber
    \end{align}

where the :math:`q_j` are linear combinations of external momenta :math:`p_i` and loop momenta :math:`k_l`.

In the first step of our approach, the loop integral is 
converted from the momentum representation to the Feynman parameter representation, see for example [Hei08]_ (Chapter 3).

The module :mod:`pySecDec.loop_integral` implements exactly that conversion.
The most basic use is to calculate the first and the second 
Symanzik polynomial ``U`` and ``F``, respectively, from the propagators of a loop integral.

.. _one-loop-bubble:

One Loop Bubble
~~~~~~~~~~~~~~~

To calculate ``U`` and ``F`` of the one loop bubble, type the following
commands:

>>> from pySecDec.loop_integral import LoopIntegralFromPropagators
>>> propagators = ['k**2', '(k - p)**2']
>>> loop_momenta = ['k']
>>> one_loop_bubble = LoopIntegralFromPropagators(propagators, loop_momenta)
>>> one_loop_bubble.U
 + (1)*x0 + (1)*x1
>>> one_loop_bubble.F
 + (-p**2)*x0*x1

The example above among other useful features is also stated in the full
documentation of
:class:`LoopIntegralFromPropagators() <pySecDec.loop_integral.LoopIntegralFromPropagators>`
in the reference guide.

Two Loop Planar Box with Numerator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Consider the propagators of the two loop planar box:

>>> propagators = ['k1**2','(k1+p2)**2',
...                '(k1-p1)**2','(k1-k2)**2',
...                '(k2+p2)**2','(k2-p1)**2',
...                '(k2+p2+p3)**2']
>>> loop_momenta = ['k1','k2']

We could now instantiate the :class:`LoopIntegral <pySecDec.loop_integral.LoopIntegral>`
just like :ref:`before <one-loop-bubble>`. However, let us consider an additional numerator:

>>> numerator = 'k1(mu)*k1(mu) + 2*k1(mu)*p3(mu) + p3(mu)*p3(mu)' # (k1 + p3) ** 2

In order to unambiguously define the loop integral, we must state which
symbols denote the ``Lorentz_indices``
(just ``mu`` in this case here) and the ``external_momenta``:

>>> external_momenta = ['p1','p2','p3','p4']
>>> Lorentz_indices=['mu']

With that, we can Feynman parametrize the two loop box with a numerator:

>>> box = LoopIntegralFromPropagators(propagators, loop_momenta, external_momenta,
...                                     numerator=numerator, Lorentz_indices=Lorentz_indices)
>>> box.U
 + (1)*x3*x6 + (1)*x3*x5 + (1)*x3*x4 + (1)*x2*x6 + (1)*x2*x5 + (1)*x2*x4 + (1)*x2*x3 + (1)*x1*x6 + (1)*x1*x5 + (1)*x1*x4 + (1)*x1*x3 + (1)*x0*x6 + (1)*x0*x5 + (1)*x0*x4 + (1)*x0*x3
>>> box.F
 + (-p1**2 - 2*p1*p2 - 2*p1*p3 - p2**2 - 2*p2*p3 - p3**2)*x3*x5*x6 + (-p3**2)*x3*x4*x6 + (-p1**2 - 2*p1*p2 - p2**2)*x3*x4*x5 + (-p1**2 - 2*p1*p2 - 2*p1*p3 - p2**2 - 2*p2*p3 - p3**2)*x2*x5*x6 + (-p3**2)*x2*x4*x6 + (-p1**2 - 2*p1*p2 - p2**2)*x2*x4*x5 + (-p1**2 - 2*p1*p2 - 2*p1*p3 - p2**2 - 2*p2*p3 - p3**2)*x2*x3*x6 + (-p1**2 - 2*p1*p2 - p2**2)*x2*x3*x4 + (-p1**2 - 2*p1*p2 - 2*p1*p3 - p2**2 - 2*p2*p3 - p3**2)*x1*x5*x6 + (-p3**2)*x1*x4*x6 + (-p1**2 - 2*p1*p2 - p2**2)*x1*x4*x5 + (-p3**2)*x1*x3*x6 + (-p1**2 - 2*p1*p2 - p2**2)*x1*x3*x5 + (-p1**2 - 2*p1*p2 - p2**2)*x1*x2*x6 + (-p1**2 - 2*p1*p2 - p2**2)*x1*x2*x5 + (-p1**2 - 2*p1*p2 - p2**2)*x1*x2*x4 + (-p1**2 - 2*p1*p2 - p2**2)*x1*x2*x3 + (-p1**2 - 2*p1*p2 - 2*p1*p3 - p2**2 - 2*p2*p3 - p3**2)*x0*x5*x6 + (-p3**2)*x0*x4*x6 + (-p1**2 - 2*p1*p2 - p2**2)*x0*x4*x5 + (-p2**2 - 2*p2*p3 - p3**2)*x0*x3*x6 + (-p1**2)*x0*x3*x5 + (-p2**2)*x0*x3*x4 + (-p1**2)*x0*x2*x6 + (-p1**2)*x0*x2*x5 + (-p1**2)*x0*x2*x4 + (-p1**2)*x0*x2*x3 + (-p2**2)*x0*x1*x6 + (-p2**2)*x0*x1*x5 + (-p2**2)*x0*x1*x4 + (-p2**2)*x0*x1*x3
>>> box.numerator
 + (2*eps*p3(mu)**2 + 2*p3(mu)**2)*U**2 + (eps - 2)*x6*F + (eps - 2)*x5*F + (eps - 2)*x4*F + (eps - 2)*x3*F + (-4*eps*p2(mu)*p3(mu) - 4*eps*p3(mu)**2 - 4*p2(mu)*p3(mu) - 4*p3(mu)**2)*x3*x6*U + (4*eps*p1(mu)*p3(mu) + 4*p1(mu)*p3(mu))*x3*x5*U + (-4*eps*p2(mu)*p3(mu) - 4*p2(mu)*p3(mu))*x3*x4*U + (2*eps*p2(mu)**2 + 4*eps*p2(mu)*p3(mu) + 2*eps*p3(mu)**2 + 2*p2(mu)**2 + 4*p2(mu)*p3(mu) + 2*p3(mu)**2)*x3**2*x6**2 + (-4*eps*p1(mu)*p2(mu) - 4*eps*p1(mu)*p3(mu) - 4*p1(mu)*p2(mu) - 4*p1(mu)*p3(mu))*x3**2*x5*x6 + (2*eps*p1(mu)**2 + 2*p1(mu)**2)*x3**2*x5**2 + (4*eps*p2(mu)**2 + 4*eps*p2(mu)*p3(mu) + 4*p2(mu)**2 + 4*p2(mu)*p3(mu))*x3**2*x4*x6 + (-4*eps*p1(mu)*p2(mu) - 4*p1(mu)*p2(mu))*x3**2*x4*x5 + (2*eps*p2(mu)**2 + 2*p2(mu)**2)*x3**2*x4**2 + (4*eps*p1(mu)*p3(mu) + 4*p1(mu)*p3(mu))*x2*x6*U + (4*eps*p1(mu)*p3(mu) + 4*p1(mu)*p3(mu))*x2*x5*U + (4*eps*p1(mu)*p3(mu) + 4*p1(mu)*p3(mu))*x2*x4*U + (4*eps*p1(mu)*p3(mu) + 4*p1(mu)*p3(mu))*x2*x3*U + (-4*eps*p1(mu)*p2(mu) - 4*eps*p1(mu)*p3(mu) - 4*p1(mu)*p2(mu) - 4*p1(mu)*p3(mu))*x2*x3*x6**2 + (4*eps*p1(mu)**2 - 4*eps*p1(mu)*p2(mu) - 4*eps*p1(mu)*p3(mu) + 4*p1(mu)**2 - 4*p1(mu)*p2(mu) - 4*p1(mu)*p3(mu))*x2*x3*x5*x6 + (4*eps*p1(mu)**2 + 4*p1(mu)**2)*x2*x3*x5**2 + (-8*eps*p1(mu)*p2(mu) - 4*eps*p1(mu)*p3(mu) - 8*p1(mu)*p2(mu) - 4*p1(mu)*p3(mu))*x2*x3*x4*x6 + (4*eps*p1(mu)**2 - 4*eps*p1(mu)*p2(mu) + 4*p1(mu)**2 - 4*p1(mu)*p2(mu))*x2*x3*x4*x5 + (-4*eps*p1(mu)*p2(mu) - 4*p1(mu)*p2(mu))*x2*x3*x4**2 + (-4*eps*p1(mu)*p2(mu) - 4*eps*p1(mu)*p3(mu) - 4*p1(mu)*p2(mu) - 4*p1(mu)*p3(mu))*x2*x3**2*x6 + (4*eps*p1(mu)**2 + 4*p1(mu)**2)*x2*x3**2*x5 + (-4*eps*p1(mu)*p2(mu) - 4*p1(mu)*p2(mu))*x2*x3**2*x4 + (2*eps*p1(mu)**2 + 2*p1(mu)**2)*x2**2*x6**2 + (4*eps*p1(mu)**2 + 4*p1(mu)**2)*x2**2*x5*x6 + (2*eps*p1(mu)**2 + 2*p1(mu)**2)*x2**2*x5**2 + (4*eps*p1(mu)**2 + 4*p1(mu)**2)*x2**2*x4*x6 + (4*eps*p1(mu)**2 + 4*p1(mu)**2)*x2**2*x4*x5 + (2*eps*p1(mu)**2 + 2*p1(mu)**2)*x2**2*x4**2 + (4*eps*p1(mu)**2 + 4*p1(mu)**2)*x2**2*x3*x6 + (4*eps*p1(mu)**2 + 4*p1(mu)**2)*x2**2*x3*x5 + (4*eps*p1(mu)**2 + 4*p1(mu)**2)*x2**2*x3*x4 + (2*eps*p1(mu)**2 + 2*p1(mu)**2)*x2**2*x3**2 + (-4*eps*p2(mu)*p3(mu) - 4*p2(mu)*p3(mu))*x1*x6*U + (-4*eps*p2(mu)*p3(mu) - 4*p2(mu)*p3(mu))*x1*x5*U + (-4*eps*p2(mu)*p3(mu) - 4*p2(mu)*p3(mu))*x1*x4*U + (-4*eps*p2(mu)*p3(mu) - 4*p2(mu)*p3(mu))*x1*x3*U + (4*eps*p2(mu)**2 + 4*eps*p2(mu)*p3(mu) + 4*p2(mu)**2 + 4*p2(mu)*p3(mu))*x1*x3*x6**2 + (-4*eps*p1(mu)*p2(mu) + 4*eps*p2(mu)**2 + 4*eps*p2(mu)*p3(mu) - 4*p1(mu)*p2(mu) + 4*p2(mu)**2 + 4*p2(mu)*p3(mu))*x1*x3*x5*x6 + (-4*eps*p1(mu)*p2(mu) - 4*p1(mu)*p2(mu))*x1*x3*x5**2 + (8*eps*p2(mu)**2 + 4*eps*p2(mu)*p3(mu) + 8*p2(mu)**2 + 4*p2(mu)*p3(mu))*x1*x3*x4*x6 + (-4*eps*p1(mu)*p2(mu) + 4*eps*p2(mu)**2 - 4*p1(mu)*p2(mu) + 4*p2(mu)**2)*x1*x3*x4*x5 + (4*eps*p2(mu)**2 + 4*p2(mu)**2)*x1*x3*x4**2 + (4*eps*p2(mu)**2 + 4*eps*p2(mu)*p3(mu) + 4*p2(mu)**2 + 4*p2(mu)*p3(mu))*x1*x3**2*x6 + (-4*eps*p1(mu)*p2(mu) - 4*p1(mu)*p2(mu))*x1*x3**2*x5 + (4*eps*p2(mu)**2 + 4*p2(mu)**2)*x1*x3**2*x4 + (-4*eps*p1(mu)*p2(mu) - 4*p1(mu)*p2(mu))*x1*x2*x6**2 + (-8*eps*p1(mu)*p2(mu) - 8*p1(mu)*p2(mu))*x1*x2*x5*x6 + (-4*eps*p1(mu)*p2(mu) - 4*p1(mu)*p2(mu))*x1*x2*x5**2 + (-8*eps*p1(mu)*p2(mu) - 8*p1(mu)*p2(mu))*x1*x2*x4*x6 + (-8*eps*p1(mu)*p2(mu) - 8*p1(mu)*p2(mu))*x1*x2*x4*x5 + (-4*eps*p1(mu)*p2(mu) - 4*p1(mu)*p2(mu))*x1*x2*x4**2 + (-8*eps*p1(mu)*p2(mu) - 8*p1(mu)*p2(mu))*x1*x2*x3*x6 + (-8*eps*p1(mu)*p2(mu) - 8*p1(mu)*p2(mu))*x1*x2*x3*x5 + (-8*eps*p1(mu)*p2(mu) - 8*p1(mu)*p2(mu))*x1*x2*x3*x4 + (-4*eps*p1(mu)*p2(mu) - 4*p1(mu)*p2(mu))*x1*x2*x3**2 + (2*eps*p2(mu)**2 + 2*p2(mu)**2)*x1**2*x6**2 + (4*eps*p2(mu)**2 + 4*p2(mu)**2)*x1**2*x5*x6 + (2*eps*p2(mu)**2 + 2*p2(mu)**2)*x1**2*x5**2 + (4*eps*p2(mu)**2 + 4*p2(mu)**2)*x1**2*x4*x6 + (4*eps*p2(mu)**2 + 4*p2(mu)**2)*x1**2*x4*x5 + (2*eps*p2(mu)**2 + 2*p2(mu)**2)*x1**2*x4**2 + (4*eps*p2(mu)**2 + 4*p2(mu)**2)*x1**2*x3*x6 + (4*eps*p2(mu)**2 + 4*p2(mu)**2)*x1**2*x3*x5 + (4*eps*p2(mu)**2 + 4*p2(mu)**2)*x1**2*x3*x4 + (2*eps*p2(mu)**2 + 2*p2(mu)**2)*x1**2*x3**2

We can also generate the output in terms of Mandelstam invariants:

>>> replacement_rules = [
...                        ('p1*p1', 0),
...                        ('p2*p2', 0),
...                        ('p3*p3', 0),
...                        ('p4*p4', 0),
...                        ('p1*p2', 's/2'),
...                        ('p2*p3', 't/2'),
...                        ('p1*p3', '-s/2-t/2')
...                     ]
>>> box = LoopIntegralFromPropagators(propagators, loop_momenta, external_momenta,
...                                     numerator=numerator, Lorentz_indices=Lorentz_indices,
...                                     replacement_rules=replacement_rules)
>>> box.U
 + (1)*x3*x6 + (1)*x3*x5 + (1)*x3*x4 + (1)*x2*x6 + (1)*x2*x5 + (1)*x2*x4 + (1)*x2*x3 + (1)*x1*x6 + (1)*x1*x5 + (1)*x1*x4 + (1)*x1*x3 + (1)*x0*x6 + (1)*x0*x5 + (1)*x0*x4 + (1)*x0*x3
>>> box.F
 + (-s)*x3*x4*x5 + (-s)*x2*x4*x5 + (-s)*x2*x3*x4 + (-s)*x1*x4*x5 + (-s)*x1*x3*x5 + (-s)*x1*x2*x6 + (-s)*x1*x2*x5 + (-s)*x1*x2*x4 + (-s)*x1*x2*x3 + (-s)*x0*x4*x5 + (-t)*x0*x3*x6
>>> box.numerator
 + (eps - 2)*x6*F + (eps - 2)*x5*F + (eps - 2)*x4*F + (eps - 2)*x3*F + (-2*eps*t - 2*t)*x3*x6*U + (-2*eps*s - 2*eps*t - 2*s - 2*t)*x3*x5*U + (-2*eps*t - 2*t)*x3*x4*U + (2*eps*t + 2*t)*x3**2*x6**2 + (2*eps*t + 2*t)*x3**2*x5*x6 + (2*eps*t + 2*t)*x3**2*x4*x6 + (-2*eps*s - 2*s)*x3**2*x4*x5 + (-2*eps*s - 2*eps*t - 2*s - 2*t)*x2*x6*U + (-2*eps*s - 2*eps*t - 2*s - 2*t)*x2*x5*U + (-2*eps*s - 2*eps*t - 2*s - 2*t)*x2*x4*U + (-2*eps*s - 2*eps*t - 2*s - 2*t)*x2*x3*U + (2*eps*t + 2*t)*x2*x3*x6**2 + (2*eps*t + 2*t)*x2*x3*x5*x6 + (-2*eps*s + 2*eps*t - 2*s + 2*t)*x2*x3*x4*x6 + (-2*eps*s - 2*s)*x2*x3*x4*x5 + (-2*eps*s - 2*s)*x2*x3*x4**2 + (2*eps*t + 2*t)*x2*x3**2*x6 + (-2*eps*s - 2*s)*x2*x3**2*x4 + (-2*eps*t - 2*t)*x1*x6*U + (-2*eps*t - 2*t)*x1*x5*U + (-2*eps*t - 2*t)*x1*x4*U + (-2*eps*t - 2*t)*x1*x3*U + (2*eps*t + 2*t)*x1*x3*x6**2 + (-2*eps*s + 2*eps*t - 2*s + 2*t)*x1*x3*x5*x6 + (-2*eps*s - 2*s)*x1*x3*x5**2 + (2*eps*t + 2*t)*x1*x3*x4*x6 + (-2*eps*s - 2*s)*x1*x3*x4*x5 + (2*eps*t + 2*t)*x1*x3**2*x6 + (-2*eps*s - 2*s)*x1*x3**2*x5 + (-2*eps*s - 2*s)*x1*x2*x6**2 + (-4*eps*s - 4*s)*x1*x2*x5*x6 + (-2*eps*s - 2*s)*x1*x2*x5**2 + (-4*eps*s - 4*s)*x1*x2*x4*x6 + (-4*eps*s - 4*s)*x1*x2*x4*x5 + (-2*eps*s - 2*s)*x1*x2*x4**2 + (-4*eps*s - 4*s)*x1*x2*x3*x6 + (-4*eps*s - 4*s)*x1*x2*x3*x5 + (-4*eps*s - 4*s)*x1*x2*x3*x4 + (-2*eps*s - 2*s)*x1*x2*x3**2

.. _sector_decomposition:

Sector Decomposition
--------------------

The sector decomposition algorithm aims to factorize the polynomials :math:`P_i`
as products of a monomial and a polynomial with nonzero constant term:

.. math::

    P_i( \{x_j\} ) \longmapsto \prod_j x_j^{\alpha_j} \left( const + p_i(\{ x_j \}) \right).

Factorizing polynomials in that way by expoliting integral transformations
is the first step in an algorithm for solving dimensionally
regulated integrals of the form

.. math::

    \int_0^1 \prod_{i,j} P_i(\{ x_j \})^{\beta_i} ~ dx_j.

The iterative sector decomposition splits the integral and remaps the integration domain
until all polynomials :math:`P_i` in all arising integrals (called `sectors`) have the
desired form :math:`const + polynomial`.
An introduction to the sector decomposition approach can be found in [Hei08]_.

To demonstrate the :mod:`pySecDec.decomposition` module, we decompose the polynomials

>>> p1 = Polynomial.from_expression('x + A*y', ['x','y','z'])
>>> p2 = Polynomial.from_expression('x + B*y*z', ['x','y','z'])

Let us first focus on the iterative decomposition of ``p1``. In the `pySecDec`
framework, we first have to pack ``p1`` into a :class:`Sector <pySecDec.decomposition.Sector>`:

>>> from pySecDec.decomposition import Sector
>>> initial_sector = Sector([p1])
>>> print(initial_sector)
Sector:
Jacobian= + (1)
cast=[( + (1)) * ( + (1)*x + (A)*y)]
other=[]

We can now run the iterative decomposition and take a look at the decomposed
sectors:

.. code:: python

    >>> from pySecDec.decomposition.iterative import iterative_decomposition
    >>> decomposed_sectors = iterative_decomposition(initial_sector)
    >>> for sector in decomposed_sectors:
    ...     print(sector)
    ...     print('\n')
    ...
    Sector:
    Jacobian= + (1)*x
    cast=[( + (1)*x) * ( + (1) + (A)*y)]
    other=[]


    Sector:
    Jacobian= + (1)*y
    cast=[( + (1)*y) * ( + (1)*x + (A))]
    other=[]


The decomposition of ``p2`` needs two iterations and yields three sectors:

.. code:: python

    >>> initial_sector = Sector([p2])
    >>> decomposed_sectors = iterative_decomposition(initial_sector)
    >>> for sector in decomposed_sectors:
    ...     print(sector)
    ...     print('\n')
    ...
    Sector:
    Jacobian= + (1)*x
    cast=[( + (1)*x) * ( + (1) + (B)*y*z)]
    other=[]


    Sector:
    Jacobian= + (1)*x*y
    cast=[( + (1)*x*y) * ( + (1) + (B)*z)]
    other=[]


    Sector:
    Jacobian= + (1)*y*z
    cast=[( + (1)*y*z) * ( + (1)*x + (B))]
    other=[]


Note that we declared ``z`` as a variable for sector ``p1`` evne though it does not depend on it.
This declaration is necessary if we want to simultaneously decompose ``p1`` and ``p2``:


.. code:: python

    >>> initial_sector = Sector([p1, p2])
    >>> decomposed_sectors = iterative_decomposition(initial_sector)
    >>> for sector in decomposed_sectors:
    ...      print(sector)
    ...      print('\n')
    ...
    Sector:
    Jacobian= + (1)*x
    cast=[( + (1)*x) * ( + (1) + (A)*y), ( + (1)*x) * ( + (1) + (B)*y*z)]
    other=[]


    Sector:
    Jacobian= + (1)*x*y
    cast=[( + (1)*y) * ( + (1)*x + (A)), ( + (1)*x*y) * ( + (1) + (B)*z)]
    other=[]


    Sector:
    Jacobian= + (1)*y*z
    cast=[( + (1)*y) * ( + (1)*x*z + (A)), ( + (1)*y*z) * ( + (1)*x + (B))]
    other=[]


We just fully decomposed ``p1`` and ``p2``. In some cases, one may want to bring
one polynomial, say ``p1``, into standard form, but not necessarily the other.
For that purpose, the :class:`Sector <pySecDec.decomposition.Sector>` can take
a second argument. In the following code example, we bring ``p1`` into standard
form, apply all transformations to ``p2`` as well, but stop before ``p2`` is fully
decomposed:


.. code:: python

    >>> initial_sector = Sector([p1], [p2])
    >>> decomposed_sectors = iterative_decomposition(initial_sector)
    >>> for sector in decomposed_sectors:
    ...      print(sector)
    ...      print('\n')
    ...
    Sector:
    Jacobian= + (1)*x
    cast=[( + (1)*x) * ( + (1) + (A)*y)]
    other=[ + (1)*x + (B)*x*y*z]


    Sector:
    Jacobian= + (1)*y
    cast=[( + (1)*y) * ( + (1)*x + (A))]
    other=[ + (1)*x*y + (B)*y*z]


Subtraction
-----------

In the subtraction, we want to perform those integrations
that lead to :math:`\epsilon` divergencies. The master formula
for one integration variables is

.. math::
    \int_0^1 {x^{(a - b \epsilon)} \mathcal{I} (x, \epsilon) dx} =
    \sum_{p=0}^{|a|-1}
        {
            \frac{1}{a + p + 1 - b \epsilon}
            \frac{\mathcal{I}^{(p)} (0, \epsilon)}{p!} +
            \int_0^1
            {
                x^{(a - b \epsilon)}
                R(x, \epsilon) dx
            }
        }

where :math:`\mathcal{I}^{(p)}` is denotes the p-th derivative
of :math:`\mathcal{I}` with respect to :math:`x`. The equation
above effectively defines the remainder term :math:`R`.
All terms on the right hand side of the equation above are
constructed to be free of divergencies. For more details
and the generalization to multiple variables, we refer the
reader to [Hei08]_.
In the following, we show how to use the implementation in
`pySecDec`.

To initialize the subtraction, we first define a factorized
expression of the form
:math:`x^{(-1 - b_x \epsilon)} y^{(-2 - b_y \epsilon)} \mathcal{I} (x, y, \epsilon)`:

>>> from pySecDec.algebra import Expression
>>> symbols = ['x','y','eps']
>>> x_monomial = Expression('x**(-1 - b_x*eps)', symbols)
>>> y_monomial = Expression('y**(-2 - b_y*eps)', symbols)
>>> cal_I = Expression('cal_I(x, y, eps)', symbols)

We must pack the monomials into a :class:`pySecDec.algebra.Product`:

>>> from pySecDec.algebra import Product
>>> monomials = Product(x_monomial, y_monomial)

Although this seems to be to complete input according to the equation
above, we are still missing a structure to store poles in. The function
:func:`pySecDec.subtraction.integrate_pole_part` is designed to return
an iterable of the same type as the input. That is particularly important
since the output of the subtraction of one variable is the input for the
subtraction of the next variable. We will see this iteration later. Initially,
we do not have poles yet, therefore we define a `one` of the required type:

>>> from pySecDec.algebra import Pow
>>> import numpy as np
>>> polynomial_one = Polynomial(np.zeros([1,len(symbols)], dtype=int), np.array([1]), symbols, copy=False)
>>> pole_part_initializer = Pow(polynomial_one, -polynomial_one)

``pole_part_initializer`` is of type :class:`pySecDec.algebra.Pow` and has ``-polynomial_one``
in the exponent. We initialize the `base` with ``polynomial_one``; i.e. a one packed into
a polynomial. The function :func:`pySecDec.subtraction.integrate_pole_part` populates the
`base` with factors of :math:`b\epsilon` when poles arise.

We are now ready to build the ``subtraction_initializer`` - the :class:`pySecDec.algebra.Product`
to be passed into :func:`pySecDec.subtraction.integrate_pole_part`.

>>> from pySecDec.subtraction import integrate_pole_part
>>> subtraction_initializer = Product(monomials, pole_part_initializer, cal_I)
>>> x_subtracted = integrate_pole_part(subtraction_initializer, 0)

The second argument of :func:`pySecDec.subtraction.integrate_pole_part` specifies
to which variable we want to apply the master formula, here we choose :math:`x`.
First, remember that the x monomial is a dimensionally regulated :math:`x^-1`.
Therefore, the sum collapses to only one term and we have two terms in total.
Each term corresponds to one entry in the list ``x_subtracted``:

>>> len(x_subtracted)
2

``x_subtracted`` has the same structure as our input. The first factor of each term
stores the remaining monomials:

>>> x_subtracted[0].factors[0]
(( + (1))**( + (-b_x)*eps + (-1))) * (( + (1)*y)**( + (-b_y)*eps + (-2)))
>>> x_subtracted[1].factors[0]
(( + (1)*x)**( + (-b_x)*eps + (-1))) * (( + (1)*y)**( + (-b_y)*eps + (-2)))

The second factor stores the :math:`\epsilon` poles. There is an epsilon pole in the first term, but
still none in the second:

>>> x_subtracted[0].factors[1]
( + (-b_x)*eps) ** ( + (-1))
>>> x_subtracted[1].factors[1]
( + (1)) ** ( + (-1))

The last factor catches everything that is not covered by the first two fields:

>>> x_subtracted[0].factors[2]
(cal_I( + (0), + (1)*y, + (1)*eps))
>>> x_subtracted[1].factors[2]
(cal_I( + (1)*x, + (1)*y, + (1)*eps)) + (( + (-1)) * (cal_I( + (0), + (1)*y, + (1)*eps)))

We have now performed the subtraction for :math:`x`. Because in and output have a similar
structure, we can easily perform the subtraction for :math:`y` as well:

.. code:: python

    >>> x_and_y_subtracted = []
    >>> for s in x_subtracted:
    ...     x_and_y_subtracted.extend( integrate_pole_part(s,1) )

Alternatively, we can directly instruct :func:`pySecDec.subtraction.integrate_pole_part`
to perform both subtractions:

>>> alternative_x_and_y_subtracted = integrate_pole_part(subtraction_initializer,0,1)

In both cases, the result is a list of the terms appearing on the right hand side of the
master equation.

Expansion
---------

The purpose of the :mod:`expansion <pySecDec.expansion>` module is,
as the name suggests, to provide routines to perform a series expansion.
The module basically implements two routines - the Taylor expansion
(:func:`pySecDec.expansion.expand_Taylor`) and an expansion of polyrational
functions supporting singularities in the expansion variable
(:func:`pySecDec.expansion.expand_singular`).

.. _Taylor_intro:

Taylor Expansion
~~~~~~~~~~~~~~~~

The function :func:`pySecDec.expansion.expand_Taylor` implements the ordinary
Taylor expansion. It takes an algebraic expression (in the sense of the
:ref:`algebra module <algebra_intro>`, the index of the expansion variable
and the order to which the expression shall be expanded:

>>> from pySecDec.algebra import Expression
>>> from pySecDec.expansion import expand_Taylor
>>> expression = Expression('x**eps', ['eps'])
>>> expand_Taylor(expression, 0, 2).simplify()
 + (1) + (log( + (x)))*eps + ((log( + (x))) * (log( + (x))) * ( + (1/2)))*eps**2

It is also possible to expand an expression in multiple variables simultaneously:

>>> expression = Expression('x**(eps + alpha)', ['eps', 'alpha'])
>>> expand_Taylor(expression, [0,1], [2,0]).simplify()
 + (1) + (log( + (x)))*eps + ((log( + (x))) * (log( + (x))) * ( + (1/2)))*eps**2

The command above instructs :func:`pySecDec.expansion.expand_Taylor` to expand
the ``expression`` to the second order in the variable indexed ``0`` (``eps``)
and to the zeroth order in the variable indexed ``1`` (``alpha``).

Laurent Expansion
~~~~~~~~~~~~~~~~~

:func:`pySecDec.expansion.expand_singular` Laurent expands polyrational functions.

Its input is more restrictive than for the :ref:`Taylor expansion <Taylor_intro>`.
It expects a :class:`Product <pySecDec.algebra.Product>` where the factors are either
:class:`Polynomials <pySecDec.algebra.Polynomial>` or
:class:`ExponentiatedPolynomials <pySecDec.algebra.ExponentiatedPolynomial>`
with ``exponent = -1``:

>>> from pySecDec.expansion import expand_singular
>>> expression = Expression('1/(eps + alpha)', ['eps', 'alpha']).simplify()
>>> expand_singular(expression, 0, 1)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/home/pcl340a/sjahn/Projects/pySecDec/pySecDec/expansion.py", line 241, in expand_singular
    return _expand_and_flatten(product, indices, orders, _expand_singular_step)
  File "/home/pcl340a/sjahn/Projects/pySecDec/pySecDec/expansion.py", line 209, in _expand_and_flatten
    expansion = recursive_expansion(expression, indices, orders)
  File "/home/pcl340a/sjahn/Projects/pySecDec/pySecDec/expansion.py", line 198, in recursive_expansion
    expansion = expansion_one_variable(expression, index, order)
  File "/home/pcl340a/sjahn/Projects/pySecDec/pySecDec/expansion.py", line 82, in _expand_singular_step
    raise TypeError('`product` must be a `Product`')
TypeError: `product` must be a `Product`
>>> expression # ``expression`` is indeed a polyrational function.
( + (1)*alpha + (1)*eps)**(-1)
>>> type(expression) # It is just not packed in a ``Product`` as ``expand_singular`` expects.
<class 'pySecDec.algebra.ExponentiatedPolynomial'>
>>> from pySecDec.algebra import Product
>>> expression = Product(expression)
>>> expand_singular(expression, 0, 1)
 + (( + (1)) * (( + (1)*alpha)**(-1))) + (( + (-1)) * (( + (1)*alpha**2)**(-1)))*eps

Like in the :ref:`Taylor expansion <Taylor_intro>`, we can expand simultaneously in
multiple parameters. Note, however, that the result of the Laurent expansion depends
on the ordering of the expansion variables. The second argument of :func:`pySecDec.expansion.expand_singular`
determines the order of the expansion:

>>> expression = Expression('1/(2*eps) * 1/(eps + alpha)', ['eps', 'alpha']).simplify()
>>> eps_first = expand_singular(expression, [0,1], [1,1])
>>> eps_first
 + (( + (1/2)) * (( + (1))**(-1)))*eps**-1*alpha**-1 + (( + (-1/2)) * (( + (1))**(-1)))*alpha**-2 + (( + (1)) * (( + (2))**(-1)))*eps*alpha**-3
>>> alpha_first = expand_singular(expression, [1,0], [1,1])
>>> alpha_first
 + (( + (1/2)) * (( + (1))**(-1)))*eps**-2 + (( + (-1/2)) * (( + (1))**(-1)))*eps**-3*alpha

The expression printed out by our algebra module are quite messy. In order to obtain nicer
output, we can convert these expressions to the slower but more high level `sympy`:

>>> import sympy as sp
>>> eps_first = expand_singular(expression, [0,1], [1,1])
>>> alpha_first = expand_singular(expression, [1,0], [1,1])
>>> sp.sympify(str(eps_first))
1/(2*alpha*eps) - 1/(2*alpha**2) + eps/(2*alpha**3)
>>> sp.sympify(str(alpha_first))
-alpha/(2*eps**3) + 1/(2*eps**2)
