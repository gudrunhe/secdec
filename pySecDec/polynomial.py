"""This file defines the Polynomial class"""

import numpy as np

class Polynomial(object):
    '''
    Container class for polynomials.
    Store a polynomial as list of lists counting the powers of
    the variables. For example the monomial "x1^2 + x1*x2" is
    stored as [[2,0],[1,1]].

    Coefficients are stored in a separate list of strings, e.g.
    "A*x1^2 + B*x1*x2" <-> [[2,0],[1,1]] and ["A","B"].

    :param expolist:
        iterable of iterables;
        The variable's powers for each term.

    :param coeffs:
        iterable of strings;
        The symbolic coefficients of the polynomial.

    '''
    def __init__(self, expolist, coeffs):
        self.expolist = np.array(expolist)
        assert len(self.expolist.shape) == 2, 'All entries in `expolist` must have the same length'
        if not np.issubdtype(self.expolist.dtype, np.integer):
            raise TypeError('All entries in `expolist` must be integer.')
        self.coeffs = list(coeffs)
        assert len(self.expolist) == len(self.coeffs), \
            '`expolist` (length %i) and `coeffs` (length %i) must have the same length.' %(len(self.expolist),len(self.coeffs))

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        outstr = ''
        for coeff,expolist in zip(self.coeffs,self.expolist):
            if coeff != '': outstr += " + %s" % coeff
            else:           outstr += " + 1"
            for i,power in enumerate(expolist):
                outstr += "*x%i^%i" %(i,power)

        return outstr

    def copy(self):
        "Return a copy of a :class:`.Polynomial`."
        return Polynomial(self.expolist, self.coeffs)

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
        for factor in self.factors:
            if factor.expolist.shape[1] != self.factors[0].expolist.shape[1]:
                raise TypeError('Must have the same number of variables for all factors.')

    def copy(self):
        "Return a copy of a :class:`.PolynomialProduct`."
        return PolynomialProduct(*self.factors)
