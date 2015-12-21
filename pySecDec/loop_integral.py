"""Routines to Feynman parametrize a loop integral"""

from .algebra import Polynomial
from .misc import det, adjugate
import sympy as sp
import numpy as np

class LoopIntegral(object):
    '''
    Container class for a loop integrals.
    The main purpose of this class is to convert a
    loop integral from the momentum representation
    to the Feynman parameter representation.

    It is possible to provide either the graph
    of the loop integrals as adjacency list,
    or the propagators.

    .. seealso::
        * input as graph: :meth:`.from_graph`
        * input as list of propagators: :meth:`from_propagators`

    '''
    def __init__(self, *args, **kwargs):
        if args[0] != 'using a named constructor':
            raise AttributeError('Use one of the named constructors.')

    @staticmethod
    def from_propagators(propagators, loop_momenta, external_momenta=None, symbols=None, \
                         numerator=None, replacement_rules=None, Feynman_parameters='x', \
                         regulator='eps', dimensionality='4 - 2*eps', metric='g'):
        r'''
        Construct the loop integral from a list of the
        loop momenta and a list of the propagators.

        :param propagators:
            iterable of strings or sympy expressions;
            The propagators, e.g. ['k1**2', '(k1-k2)**2 - m1**2'].

        :param loop_momenta:
            iterable of strings or sympy expressions;
            The loop momenta, e.g. ['k1','k2'].

        :param external_momenta:
            iterable of strings or sympy expressions,
            optional;
            The external momenta, e.g. ['p1','p2'].
            Specifying the `external_momenta` is only
            required when a `numerator` is to be
            constructed.

        :param symbols:
            iterable of strings or sympy expressions,
            optional;
            Any scalar symbol that appears in the `numerator`
            must be specified here, in particular all the
            Mandelstam variables and masses. If a `numerator`
            shall be constructed, all scalar products of
            external momenta that appear in the `numerator`
            must be expressed in terms scalars (the Mandelstam
            variables) in the `replacement_rules`.
            Specifying the `kinematic_invariants` is only
            required when a `numerator` is to be constructed.

        :param numerator:
            string or sympy expression, optional;
            The numerator of the loop integral.
            The numerator should be a scalar; i.e. all momenta
            should be contracted in scalar products. The scalar
            products should be in index notation e.g.
            "k1(mu)*k2(mu)". The numerator should be a sum of
            products of exclusively:
            * numbers
            * scalar products (e.g. "p1(mu)*k1(mu)*p1(nu)*k2(nu)")
            * `symbols` (e.g. "m")
            Examples:
                * 'p1(mu)*k1(mu)*p1(nu)*k2(nu) + 4*s*eps*k1(mu)*k1(mu)'
                * 'p1(mu)*(k1(mu) + k2(mu))*p1(nu)*k2(nu)'

            .. hint::
                It is possible to use numbers as indices, for example
                'p1(mu)*p2(mu)*k1(nu)*k2(nu) = p1(1)*p2(1)*k1(2)*k2(2)'.

            .. warning::
                Do **NOT** use the symbol "ScalarProduct".
                "ScalarProduct" is internally used to group the scalar
                products in index notation.

        :param replacement_rules:
            iterable of iterables with two strings or sympy
            expressions, optional;
            Symbolic replacements to be made for the external
            momenta, e.g. definition of Mandelstam variables.
            Example: [('p1*p2', 's'), ('p1**2', 0)] where
            ``p1`` and ``p2`` are external momenta.
            It is also possible to specify vector replacements,
            for example [('p4', '-(p1+p2+p3)')].

        :param Feynman_parameters:
            iterable or string, optional;
            The symbols to be used for the Feynman parameters.
            If a string is passed, the Feynman parameter
            variables will be consecutively numbered starting
            from zero.

        :param regulator:
            string or sympy symbol, optional;
            The symbol to be used for the dimensional regulator
            (typically :math:`\epsilon` or :math:`\epsilon_D`)

        :param dimensionality:
            string or sympy symbol, optional;
            The dimensionality; this also defines the symbol
            to be used for the dimensional regulator (typically
            :math:`\epsilon` or :math:`\epsilon_D`)

        :param metric:
            string or sympy symbol, optional;
            The symbol to be used for the (Minkowski) metric
            tensor :math:`g^{\mu\nu}`.

        '''
        # TODO: explain the main member properties (U, F, numerator)
        # TODO: carefully reread and check this documentation
        # TODO: implement numerator
        # TODO: check that the propagators are at most quadratic in the loop momenta
        # TODO: remove all sympy symbols from the final output U, F, and N

        self = LoopIntegral('using a named constructor')

        # sympify and store `loop_momenta`
        self.loop_momenta = []
        for expression in loop_momenta:
            expression = sp.sympify(expression)
            assert expression.is_Symbol, 'Each of the `loop_momenta` must be a `sympy symbol`.' #TODO: test error message
            self.loop_momenta.append(expression)
        self.L = len(self.loop_momenta)

        # sympify and store `propagators`
        self.propagators = sp.sympify(list(propagators))
        self.P = len(self.propagators)

        # sympify and store `Feynman_parameters`
        # There should be one Feynman parameter for each propagator.
        if isinstance(Feynman_parameters, str):
            self.Feynman_parameters = [sp.sympify(Feynman_parameters + str(i)) for i in range(self.P)]
        else:
            self.Feynman_parameters = sp.sympify(list(Feynman_parameters))
            assert len(self.Feynman_parameters) == len(self.propagators), \
                'Mismatch between the number of `propagators` (%i) and the number of `Feynman_parameters` (%i)' % \
                ( len(self.propagators) , len(self.Feynman_parameters) )

        if replacement_rules is None:
            self.replacement_rules = []
        else:
            self.replacement_rules = np.array(replacement_rules)
            assert len(self.replacement_rules.shape) == 2, "Wrong format for `replacement_rules`" #TODO: test error message
            assert self.replacement_rules.shape[1] == 2 , "Wrong format for `replacement_rules`" #TODO: test error message

        return self

    @staticmethod
    def from_graph(): #TODO: What input is required for the cut construct?
        # TODO: implement replacement rules in cut construct
        raise NotImplementedError()
        self = LoopIntegral()
        # TODO: implement cut construction
        return self

    @property
    def propagator_sum(self):
        while True:
            try:
                return self._propagator_sum
            except AttributeError:
                self._propagator_sum = sum(prop*FP for prop,FP in zip(self.propagators,self.Feynman_parameters)).expand()

    @property
    def M(self):
        while True:
            try:
                return self._M
            except AttributeError:
                M = np.empty((self.L,self.L), dtype=object)
                for i in range(self.L):
                    for j in range(i,self.L):
                        current_term = self.propagator_sum.coeff(self.loop_momenta[i]*self.loop_momenta[j])
                        if i != j:
                            current_term /= 2
                        tmp = Polynomial.from_expression(current_term, self.Feynman_parameters)
                        # all coeffs of the polynomials in M must be integer
                        # convert to `int` since it is faster than calculating with sympy expressions
                        tmp.coeffs = tmp.coeffs.astype(int)
                        # M is symmetric
                        M[j,i] = M[i,j] = tmp
                self._M = M

    @property
    def detM(self):
        while True:
            try:
                return self._detM
            except AttributeError:
                self._detM = det(self.M)

    @property
    def aM(self):
        while True:
            try:
                return self._aM
            except AttributeError:
                self._aM = adjugate(self.M)

    @property
    def Q(self):
        while True:
            try:
                return self._Q
            except AttributeError:
                Q = np.empty(self.L, dtype=object)
                for i in range(self.L):
                    current_term = self.propagator_sum.coeff(self.loop_momenta[i]).subs([(l,0) for l in self.loop_momenta])/(-2)
                    Q[i] = Polynomial.from_expression(current_term.expand().subs(self.replacement_rules), self.Feynman_parameters)
                self._Q = Q

    @property
    def J(self):
        while True:
            try:
                return self._J
            except AttributeError:
                sympy_J = self.propagator_sum.subs([(l,0) for l in self.loop_momenta])
                self._J = Polynomial.from_expression(sympy_J.expand().subs(self.replacement_rules), self.Feynman_parameters)

    # equation (8) of arXiv:0803.4177: U = det(M)
    U = detM

    @property
    def F(self):
        while True:
            try:
                return self._F
            except AttributeError:
                # equation (8) of arXiv:0803.4177: F = det(M)*(Q.transpose*inverse(M)*Q-J) = (Q.transpose*adjugate(M)*Q-U*J)
                F = 0
                for i in range(self.L):
                    for j in range(self.L):
                        F += self.Q[i]*self.aM[i,j]*self.Q[j]
                self._F = F - self.U * self.J
                for i,coeff in enumerate(self._F.coeffs):
                    self._F.coeffs[i] = coeff.expand().subs(self.replacement_rules)
                self._F = self._F.simplify()
