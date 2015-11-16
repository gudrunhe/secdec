"""This file defines the Polynomial class"""

import numpy as np
import sympy as sp

class Polynomial(object):
    '''
    Container class for polynomials.
    Store a polynomial as list of lists counting the powers of
    the variables. For example the polynomial "x1**2 + x1*x2" is
    stored as [[2,0],[1,1]].

    Coefficients are stored in a separate list of strings, e.g.
    "A*x0**2 + B*x0*x1" <-> [[2,0],[1,1]] and ["A","B"].

    :param expolist:
        iterable of iterables;
        The variable's powers for each term.

    :param coeffs:
        iterable;
        The coefficients of the polynomial.

    :param polysymbols:
        iterable or string, optional;
        The symbols to be used for the polynomial variables
        when converted to string. If a string is passed, the
        variables will be consecutively numbered.

        For example: expolist=[[2,0],[1,1]] coeffs=["A","B"]
         * polysymbols='x' (default) <-> "A*x0**2 + B*x0*x1"
         * polysymbols=['x','y']     <-> "A*x**2 + B*x*y"

    '''
    def __init__(self, expolist, coeffs, polysymbols='x'):
        self.expolist = np.array(expolist)
        assert len(self.expolist.shape) == 2, 'All entries in `expolist` must have the same length'
        if not np.issubdtype(self.expolist.dtype, np.integer):
            raise TypeError('All entries in `expolist` must be integer.')
        self.coeffs = list(coeffs)
        assert len(self.expolist) == len(self.coeffs), \
            '`expolist` (length %i) and `coeffs` (length %i) must have the same length.' %(len(self.expolist),len(self.coeffs))
        self.number_of_variables = self.expolist.shape[1]
        if isinstance(polysymbols, str):
            self.polysymbols=[polysymbols + str(i) for i in range(self.number_of_variables)]
        else:
            self.polysymbols=list(polysymbols)

    def __repr__(self):
        outstr = ''
        for coeff,expolist in zip(self.coeffs,self.expolist):
            if coeff != '':
                outstr += (" + (%s)" % coeff)
            else:           outstr += " + (1)"
            for i,(power,symbol) in enumerate(zip(expolist,self.polysymbols)):
                if power == 0:
                    continue
                elif power == 1:
                    outstr += "*%s" % symbol
                else:
                    outstr += "*%s**%i" %(symbol,power)
        return outstr

    __str__ = __repr__

    def copy(self):
        "Return a copy of a :class:`.Polynomial` or a subclass."
        return type(self)(self.expolist, self.coeffs, self.polysymbols)

    def has_constant_term(self):
        '''
        Return True if the polynomial can be written as:

        .. math::
            const + ...

        Otherwise, return False.

        '''
        return (self.expolist == 0).all(axis=1).any()

    def becomes_zero_for(self, zero_params):
        '''
        Return True if the polynomial becomes zero if the
        parameters passed in `zero_params` are set to zero.
        Otherwise, return False.

        :param zero_params:
            iterable of integers;
            The indices of the parameters to be checked.

        '''
        return (self.expolist > 0)[:,tuple(zero_params)].any(axis=1).all()

class SNCPolynomial(Polynomial):
    '''
    "Symbolic or Numerical Coefficiented Polynomial"
    Like :class:`.Polynomial`, but with numerical or symbolic
    coefficients (`coeffs`). The coefficients are stored
    in a numpy array.

    :param expolist:
        iterable of iterables;
        The variable's powers for each term.

    :param coeffs:
        1d array-like with numerical or sympy-symbolic
        (see http://www.sympy.org/) content, e.g. [x,1,2]
        where x is a sympy symbol;
        The coefficients of the polynomial.

    :param polysymbols:
        iterable or string, optional;
        The symbols to be used for the polynomial variables
        when converted to string. If a string is passed, the
        variables will be consecutively numbered.

        For example: expolist=[[2,0],[1,1]] coeffs=["A","B"]
         * polysymbols='x' (default) <-> "A*x0**2 + B*x0*x1"
         * polysymbols=['x','y']     <-> "A*x**2 + B*x*y"

    '''
    def __init__(self, expolist, coeffs, polysymbols='x'):
        Polynomial.__init__(self, expolist, coeffs, polysymbols)
        self.coeffs = np.array(coeffs)
        assert len(self.coeffs.shape) == 1, '`coeffs` must be one-dimensional'

    @staticmethod
    def from_expression(expression, polysymbols):
        '''
        Alternative constructor.
        Construct the polynomial from an algebraic expression.

        :param expression:
            string or sympy expression;
            The algebraic representation of the polynomial, e.g.
            "5*x1**2 + x1*x2"

        :param polysymbols:
            iterable of strings or sympy symbols;
            The symbols to be interpreted as the polynomial variables,
            e.g. "['x1','x2']".

        '''
        polysymbols = list(polysymbols)

        if not polysymbols:
            raise TypeError("`polysymbols` must contain at least one symbol")

        expression, polysymbols = sp.sympify((expression, polysymbols))

        for symbol in polysymbols:
            if not symbol.is_Symbol:
                raise TypeError("'%s' is not a symbol" % symbol)

        sympy_poly = sp.poly(expression, polysymbols)
        expolist = sympy_poly.monoms()
        coeffs = sympy_poly.coeffs()
        return SNCPolynomial(expolist, coeffs, polysymbols)

    def __add__(self, other):
        'addition operator'
        return self._sub_or_add(other, False)

    def __sub__(self, other):
        'subtraction operator'
        return self._sub_or_add(other, True)

    def _sub_or_add(self, other, sub):
        '''
        Function that implements addition and subtraction.
        The coefficients of `other` are negated if `sub`
        is `True`.

        '''
        if  type(other) is not SNCPolynomial:
            return NotImplemented

        assert self.number_of_variables == other.number_of_variables, 'Number of varibales must be equal for both polynomials in +'

        sum_expolist = np.vstack([self.expolist, other.expolist])
        sum_coeffs = np.hstack([self.coeffs, -other.coeffs if sub else other.coeffs])

        result = SNCPolynomial(sum_expolist, sum_coeffs, self.polysymbols)
        result.combine()
        return result

    def __mul__(self, other):
        'multiplication operator'
        if  type(other) is not SNCPolynomial:
            return NotImplemented

        assert self.number_of_variables == other.number_of_variables, 'Number of varibales must be equal for both factors in *'

        product_expolist = np.vstack([other.expolist + term for term in self.expolist])
        product_coeffs = np.hstack([other.coeffs * term for term in self.coeffs])

        result = SNCPolynomial(product_expolist, product_coeffs, self.polysymbols)
        result.combine()
        return result

    def __neg__(self):
        'arithmetic negation "-self"'
        return SNCPolynomial(self.expolist, [-coeff for coeff in self.coeffs], self.polysymbols)

    def combine(self):
        '''
        Combine terms that have the same exponents of
        the variables.

        '''
        for i in range(len(self.coeffs)):
            # do not have to consider terms with zero coefficient
            if self.coeffs[i] == 0: continue
            # search `self.expolist` for the same term
            same_exponents = np.where( (self.expolist[i+1:] == self.expolist[i]).all(axis=1) )
            # add all these coefficients together
            self.coeffs[i] += np.sum(self.coeffs[i+1:][same_exponents])
            # mark other terms for removal by setting coefficients to zero
            self.coeffs[i+1:][same_exponents] = 0

        # remove terms with zero coefficient
        zero_coeffs = np.where(self.coeffs == 0)
        self.coeffs = np.delete(self.coeffs, zero_coeffs)
        self.expolist = np.delete(self.expolist, zero_coeffs ,axis=0)

        # need to have at least one term
        if len(self.coeffs) == 0:
            self.coeffs = np.array([0])
            self.expolist = np.array([[0]*self.number_of_variables])

class ExponentiatedPolynomial(Polynomial):
    '''
    Like :class:`.Polynomial`, but with a global exponent.
    :math:`polynomial^{exponent}`

    :param expolist:
        iterable of iterables;
        The variable's powers for each term.

    :param coeffs:
        iterable;
        The coefficients of the polynomial.

    :param exponent:
        object, optional;
        The global exponent.

    :param polysymbols:
        iterable or string, optional;
        The symbols to be used for the polynomial variables
        when converted to string. If a string is passed, the
        variables will be consecutively numbered.

        For example: expolist=[[2,0],[1,1]] coeffs=["A","B"]
         * polysymbols='x' (default) <-> "A*x0**2 + B*x0*x1"
         * polysymbols=['x','y']     <-> "A*x**2 + B*x*y"

    '''
    def __init__(self, expolist, coeffs, exponent=1, polysymbols='x'):
        Polynomial.__init__(self, expolist, coeffs, polysymbols)
        self.exponent = exponent

    def __repr__(self):
        if self.exponent == 1:
            return super(ExponentiatedPolynomial, self).__repr__()
        else:
            return '(' + super(ExponentiatedPolynomial, self).__repr__() \
                       + ')**(%s)' % self.exponent
    __str__ = __repr__

    def copy(self):
        "Return a copy of a :class:`.Polynomial` or a subclass."
        return type(self)(self.expolist, self.coeffs, self.exponent, self.polysymbols)

class PolynomialProduct(object):
    r'''
    Product of polynomials.
    Store one or polynomials :math:`p_i` to be interpreted as
    product :math:`\prod_i p_i`.

    :param factors:
        arbitrarily many instances of :class:`.Polynomial`;
        The factors :math:`p_i`.

    :math:`p_i` can be accessed with ``self.factors[i]``.

    Example:

    .. code-block:: python

        p = PolynomialProduct(p0, p1)
        p0 = p.factors[0]
        p1 = p.factors[1]


    '''
    def __init__(self,*factors):
        self.factors = [factor.copy() for factor in factors]
        assert self.factors, 'Must have at least one factor'

        self.number_of_variables = self.factors[0].expolist.shape[1]

        for factor in self.factors:
            if factor.expolist.shape[1] != self.number_of_variables:
                raise TypeError('Must have the same number of variables for all factors.')

    def copy(self):
        "Return a copy of a :class:`.PolynomialProduct`."
        return PolynomialProduct(*self.factors)
