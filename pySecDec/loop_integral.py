"""Routines to Feynman parametrize a loop integral"""

from .algebra import Polynomial
from .misc import det, adjugate, powerset, missing, all_pairs, cached_property
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
                         numerator=1, replacement_rules=None, Feynman_parameters='x', \
                         regulator='eps', dimensionality='4-2*eps', metric_tensor='g'):
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

        :param metric_tensor:
            string or sympy symbol, optional;
            The symbol to be used for the (Minkowski) metric
            tensor :math:`g^{\mu\nu}`.

        '''
        # TODO: numerator implementation according to arXiv:1010.1667
        # TODO: allow tensor numerators --> correct documentation
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

        self.numerator_input = sp.sympify(numerator).expand()
        if self.numerator_input.is_Add:
            self.numerator_input_terms = list(self.numerator_input.args)
        else:
            self.numerator_input_terms = [self.numerator_input]

        self.metric_tensor = sp.sympify(metric_tensor) # TODO: assert this is a symbol

        self.external_momenta = sp.sympify(external_momenta) # TODO: assert all these are symbols

        return self

    @staticmethod
    def from_graph(): #TODO: What input is required for the cut construct?
        # TODO: implement replacement rules in cut construct
        raise NotImplementedError()
        self = LoopIntegral()
        # TODO: implement cut construction
        return self

    @cached_property
    def propagator_sum(self):
        return sum(prop*FP for prop,FP in zip(self.propagators,self.Feynman_parameters)).expand()

    @cached_property
    def M(self):
        M = np.empty((self.L,self.L), dtype=object)
        for i in range(self.L):
            for j in range(i,self.L):
                current_term = self.propagator_sum.coeff(self.loop_momenta[i]*self.loop_momenta[j])
                if i != j:
                    current_term /= 2
                tmp = Polynomial.from_expression(current_term, self.Feynman_parameters)
                # all coeffs of the polynomials in M must be integer
                # convert to `int` since it is faster than calculating with sympy expressions # TODO: not always possible --> try conversion but do not force it
                tmp.coeffs = tmp.coeffs.astype(int)
                # M is symmetric
                M[j,i] = M[i,j] = tmp
        return M

    @cached_property
    def detM(self):
        return det(self.M)

    @cached_property
    def aM(self):
        return adjugate(self.M)

    @cached_property
    def Q(self):
        Q = np.empty(self.L, dtype=object)
        for i in range(self.L):
            current_term = self.propagator_sum.coeff(self.loop_momenta[i]).subs([(l,0) for l in self.loop_momenta])/(-2)
            Q[i] = Polynomial.from_expression(current_term.expand().subs(self.replacement_rules), self.Feynman_parameters)
        return Q

    @cached_property
    def J(self):
        sympy_J = self.propagator_sum.subs([(l,0) for l in self.loop_momenta])
        return Polynomial.from_expression(sympy_J.expand().subs(self.replacement_rules), self.Feynman_parameters)

    # equation (8) of arXiv:0803.4177: U = det(M)
    U = detM

    @cached_property
    def F(self):
        # equation (8) of arXiv:0803.4177: F = det(M)*(Q.transpose*inverse(M)*Q-J) = (Q.transpose*adjugate(M)*Q-U*J)
        F = 0
        for i in range(self.L):
            for j in range(self.L):
                F += self.Q[i]*self.aM[i,j]*self.Q[j]
        F -= self.U * self.J
        for i,coeff in enumerate(F.coeffs):
            F.coeffs[i] = coeff.expand().subs(self.replacement_rules)
        return F.simplify()

    @cached_property
    def _numerator_tensors(self):
        '''
        Return the tensor structure of the numerator
        encoded in double indices. The return value is
        a triple of lists. In the first two lists, the
        first index of each pair is the position of the
        loop momentum in ``self.loop_momenta`` or
        ``self.external_momenta``, respectively. The second
        index is the contraction index.
        The last list contains the remaining factors of
        the corresponding term.

        Example:
            input:
                loop_momenta = ['k1', 'k2']
                external_momenta = ['p1', 'p2']
                numerator = 'k1(mu)*k2(mu) + A*k2(1)*k2(2)*p1(1)*p2(2)'

            output:
                numerator_loop_tensors = [[(0,mu),(1,mu)],[(1,1),(1,2)]]
                numerator_external_tensors = [[],[(0,1),(1,2)]]
                numerator_symbols = [1, A]
                _numerator_tensors = (numerator_loop_tensors,
                                      numerator_external_tensors,
                                      numerator_symbols)

        '''
        wildcard_index = sp.Wild('wildcard_index')
        def append_double_index(expr, lst, momenta):
            '''
            Append the double index of an expression of the form
            ``<momentum>(<index>)`` to `lst`.
            Return ``True`` if the ``<momentum>`` matches one of
            the `momenta`, otherwise return ``False``.

            '''
            for i,loop_momentum in enumerate(momenta):
                indexdict = expr.match(loop_momentum(wildcard_index))
                if indexdict is not None: # expression matches
                    lst.append((i,indexdict[wildcard_index]))
                    return True
            return False

        numerator_loop_tensors = []
        numerator_external_tensors = []
        numerator_symbols = []
        for term in self.numerator_input_terms:
            current_loop_tensor = []
            current_external_tensor = []
            current_remaining_factor = 1
            if term.is_Function: # catch special case ``term = <loop momentum>(<index>)``
                append_double_index(term, current_loop_tensor, self.loop_momenta)
            else:
                for arg in term.args:
                    # search for ``momentum(index)``
                    if arg.is_Function:
                        if not append_double_index(arg, current_loop_tensor, self.loop_momenta):
                            # continue lookup only if nothing found so far
                            if not append_double_index(arg, current_external_tensor, self.external_momenta):
                                current_remaining_factor *= arg
                    else:
                        current_remaining_factor *= arg
            numerator_loop_tensors.append(current_loop_tensor)
            numerator_external_tensors.append(current_external_tensor)
            numerator_symbols.append(current_remaining_factor)
        return numerator_loop_tensors, numerator_external_tensors, numerator_symbols

    @cached_property
    def numerator_loop_tensors(self):
        return self._numerator_tensors[0]

    @cached_property
    def numerator_external_tensors(self):
        return self._numerator_tensors[1]

    @cached_property
    def numerator_symbols(self):
        return self._numerator_tensors[2]

    @cached_property
    def numerator_index_notation(self):
        '''
        Generate the numerator in index notation according to
        section 2 in arXiv:1010.1667.

        '''
        # The implementation below conceptually follows section 2 of arXiv:1010.1667.

        # `scalar_factor` is the `r` dependent factor in the sum of equation (2.5) in arXiv:1010.1667v1
        # `scalar_factor` = `1/(-2)**(r/2)*Gamma(N_nu - dim*L/2 - r/2)*F**(r/2)
        scalar_factor = sp.symbols('scalar_factor')
        # TODO: let the user define this symbol
        # TODO: maybe use `pySecDec.algebra.Function`

        numerator = 0

        aM = self.aM
        Q = self.Q
        g = self.metric_tensor

        # `self.numerator_loop_tensors`: List of double indices for the loop momenta for each term.
        # `self.numerator_external_tensors`: List of double indices for the external momenta for each term.
        # `self.numerator_symbols`: List of other symbols for each term.
        # see also: Docstring of `self._numerator_tensors` above.
        # Loop over the terms (summands) of the numerator.
        # Each summand is expected to be a product.
        for loop_tensor, external_tensor, remainder in zip(self.numerator_loop_tensors, self.numerator_external_tensors, self.numerator_symbols):

            # Calculate the terms of the sum in equation (2.5) of arXiv:1010.1667v1:
            #  * The `scalar_factor` as explained above
            #  * The tensor "A"
            #  * The tensor "P"
            #  * The tensor of external momenta from the input (`external_tensor`)
            #  * Remaining scalar factors from the input (`remainder`)

            # Distribute the indices of the loop momenta over "A" (and "P") in all possible ways --> powerset
            # Note that "A" must carry an even number of indices --> stride=2
            for A_indices in powerset(loop_tensor, stride=2):
                r = len(A_indices)

                # "P" carries the indices "A" doesn't.
                P_indices = missing(loop_tensor, A_indices)

                # Distribute the indices of "A" over metric tensors in a completely symmetric way.
                # Note that "P" is already completely symmetric by definition.
                for tensors_A2_indices in all_pairs(A_indices):
                    # `tensors_A2_indices` corresponds to the indices in equation (2.15) in arXiv:1010.1667v1.

                    this_numerator_summand = scalar_factor(r) * remainder

                    # ---- multiply the caligraphic "A" in equation (2.5) of arXiv:1010.1667v1. ----
                    # `tensors_A2_indices` is a list of pairs of double indices describing one of the factors of equation (2.15) in arXiv:1010.1667v1.
                    # The following loop inserts the definition (2.18, arXiv:1010.1667v1) for all of these pairs.
                    for tensor_A2_indices in tensors_A2_indices:
                        this_numerator_summand *= aM[tensor_A2_indices[0][0],tensor_A2_indices[1][0]]*g(tensor_A2_indices[0][1], tensor_A2_indices[1][1])
                    # ------------------------------------------------------------------------------

                    # --------------------------- multiply the tensor "P" --------------------------
                    for external_momentum_index, Lorentz_index in P_indices:
                        this_tensor_P_factor = sp.sympify(aM[external_momentum_index].dot(Q))
                        # There should be eactly one external momentum in each term (summand) of `this_tensor_P_factor` --> attach the `Lorentz_index` to it
                        this_tensor_P_factor = this_tensor_P_factor.subs((p, p(Lorentz_index)) for p in self.external_momenta)
                        this_numerator_summand *= this_tensor_P_factor
                    # ------------------------------------------------------------------------------

                    # ----------------------- multiply the `external_tensor` -----------------------
                    for external_momentum_index, Lorentz_index in external_tensor:
                        this_numerator_summand *= self.external_momenta[external_momentum_index](Lorentz_index)
                    # ------------------------------------------------------------------------------

                    numerator += this_numerator_summand

            return numerator

            # TODO: Contract indices if possible using the `self.replacement_rules`
