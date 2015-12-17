"""Routines to Feynman parametrize a loop integral"""

from .algebra import Polynomial
from .misc import det, adjugate
import sympy as sp
import numpy as np

def uf(loop_momenta,propagators):
    r'''
    Construct the 1st (U) and 2nd (F) Symanzik Polynomials
    from a list of loop momenta and propagators.
    Return a tuple containing (U,F) as
    :class:`pySecDec.algebra.Polynomial`.

    :param loop_momenta:
       iterable of strings or sympy expression;
       The loop momenta, e.g. ['k1','k2'].

    :param propagators:
       iterable of strings or sympy expressions;
       The propagators, e.g. ['k1**2', '(k1-k2)**2 - m**2'].

    Adapted from Eq(8) arXiv:0803.4177.

    .. warning::
        Do **NOT** name any of your variables "Feynman0",
        "Feynman1", ... - these variables are internally
        reserved for the Feynman parameters.

    '''
    # convert input to sympy expressions
    loop_momenta = [sp.sympify(loop_momentum) for loop_momentum in loop_momenta]
    propagators = [sp.sympify(propagator) for propagator in propagators]

    L = len(loop_momenta)
    P = len(propagators)

    # define Feynman parameters
    Feynman_parameters = [sp.symbols('Feynman%i' % i) for i in range(P)]

    # Calculate the square bracket in equation (5) of arXiv:0803.4177v2.
    propsum = sum(prop*FP for prop,FP in zip(propagators,Feynman_parameters)).expand()


    # construct matrix M
    M = np.empty((L,L), dtype=object)
    for i in range(L):
        for j in range(i,L):
            current_term = propsum.coeff(loop_momenta[i]*loop_momenta[j])
            if i != j:
                current_term /= 2
            tmp = Polynomial.from_expression(current_term, Feynman_parameters)
            # all coeffs of the polynomials in M must be integer
            # convert to `int` since it is faster than calculating with sympy expressions
            tmp.coeffs = tmp.coeffs.astype(int)
            # M is symmetric
            M[j,i] = M[i,j] = tmp

    # construct vector Q
    Q = np.empty(L, dtype=object)
    for i in range(L):
        current_term = propsum.coeff(loop_momenta[i]).subs([(l,0) for l in loop_momenta])/(-2)
        Q[i] = Polynomial.from_expression(current_term, Feynman_parameters)

    # construct J
    sympy_J = propsum.subs([(l,0) for l in loop_momenta])
    J = Polynomial.from_expression(sympy_J, Feynman_parameters)


    # equation (8) of arXiv:0803.4177
    U = det(M)

    # generate adjugate aM of M
    aM = adjugate(M)

    # equation (8) of arXiv:0803.4177: F = det(M)*(Q.transpose*inverse(M)*Q-J) = (Q.transpose*adjugate(M)*Q-U*J)
    F = 0
    for i in range(L):
        for j in range(L):
            F += Q[i]*aM[i,j]*Q[j]
    F -= U * J

    F.polysymbols = U.polysymbols = Feynman_parameters

    return (U,F)
