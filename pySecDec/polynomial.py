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
        terms = []
        for coeff,expolist in zip(self.coeffs,self.expolist):
            if coeff != '': outstr += " + %s" % coeff
            else:           outstr += " + 1"
            for i,power in enumerate(expolist):
                outstr += "*x%i^%i" %(i,power)

        return outstr

    def copy(self):
        "Return a copy of a :class:`.Polynomial`."
        return Polynomial(self.expolist, self.coeffs)
