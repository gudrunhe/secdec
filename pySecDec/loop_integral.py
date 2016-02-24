"""Routines to Feynman parametrize a loop integral"""

from .algebra import Polynomial, ExponentiatedPolynomial
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
        raise AttributeError('Use one of the derived classes.')

    @cached_property
    def exponent_U(self):
        return len(self.propagators) - self.dimensionality / 2 * (len(self.loop_momenta) + 1) - self.highest_rank

    @property # not cached on purpose --> this is just making copies
    def exponentiated_U(self):
        return ExponentiatedPolynomial(self.U.expolist, self.U.coeffs, self.exponent_U, self.U.polysymbols)

    @cached_property
    def exponent_F(self):
        return self.dimensionality / 2 * len(self.loop_momenta) - len(self.propagators)

    @property # not cached on purpose --> this is just making copies
    def exponentiated_F(self):
        return ExponentiatedPolynomial(self.F.expolist, self.F.coeffs, self.exponent_F, self.F.polysymbols)

    @cached_property
    def Gamma_factor(self):
        # Every term factor in the sum of equation (2.5) in arXiv:1010.1667v1 comes with
        # the scalar factor `1/(-2)**(r/2)*Gamma(N_nu - dim*L/2 - r/2)*F**(r/2)`.
        # In order to keep the `numerator` free of poles in the regulator, we divide it
        # by the Gamma function with the smallest argument `N_nu - dim*L/2 - highest_rank//2`,
        # where `//` means integer division, and put it here.
        return sp.gamma(len(self.propagators) - self.dimensionality * len(self.loop_momenta)/2 - self.highest_rank//2)


class LoopIntegral_from_propagators(LoopIntegral):
    '''
        Construct the Feynman parametrization of a
        loop integral from the algebraic momentum
        representation.

        .. seealso::
            [Hei08]_, [GKR+11]_

        The Feynman parametrized integral is a product
        of the following expressions member properties:

        * ``self.regulator ** self.regulator_power``
        * ``self.Gamma_factor``
        * ``self.exponentiated_U``
        * ``self.exponentiated_F``
        * ``self.numerator``

        , where ``self`` is an instance of :class:`LoopIntegral`.

        Example:

        >>> from pySecDec.loop_integral import *
        >>> propagators = ['k**2', '(k - p)**2']
        >>> loop_momenta = ['k']
        >>> li = LoopIntegral.from_propagators(propagators, loop_momenta)
        >>> li.exponentiated_U
        ( + (1)*x0 + (1)*x1)**(2*eps - 2)
        >>> li.exponentiated_F
        ( + (-p**2)*x0*x1)**(-eps)

        The 1st (U) and 2nd (F) Symanzik polynomials and
        their exponents can also be accessed
        independently:

        >>> li.U
         + (1)*x0 + (1)*x1
        >>> li.F
         + (-p**2)*x0*x1
        >>>
        >>> li.exponent_U
        2*eps - 2
        >>> li.exponent_F
        -eps

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

        :param Lorentz_indices:
            iterable of strings or sympy expressions,
            optional;
            Symbols to be used as Lorentz indices in the
            numerator.

            .. seealso::
                parameter `numerator`

        :param numerator:
            string or sympy expression, optional;
            The numerator of the loop integral.
            Scalar products must be passed in index notation e.g.
            "k1(mu)*k2(mu)". The numerator should be a sum of
            products of exclusively:
            * numbers
            * scalar products (e.g. "p1(mu)*k1(mu)*p1(nu)*k2(nu)")
            * `symbols` (e.g. "m")

            Examples:
                * 'p1(mu)*k1(mu)*p1(nu)*k2(nu) + 4*s*eps*k1(mu)*k1(mu)'
                * 'p1(mu)*(k1(mu) + k2(mu))*p1(nu)*k2(nu)'
                * 'p1(mu)*k1(mu)*my_function(eps)'

            .. hint::
                It is possible to use numbers as indices, for example
                'p1(mu)*p2(mu)*k1(nu)*k2(nu) = p1(1)*p2(1)*k1(2)*k2(2)'.

            .. hint::
                The numerator may have uncontracted indices, e.g.
                'k1(mu)*k2(nu)'

            .. warning::
                **All** Lorentz indices (including the contracted ones)
                must be explicitly defined using the parameter
                `Lorentz_indices`.

            .. warning::
                It is assumed that the numerator is and all its
                derivatives by the `regulator` are finite and defined
                if :math:`\epsilon=0` is inserted explicitly.
                In particular, if user defined functions (like in the
                example ``p1(mu)*k1(mu)*my_function(eps)``) appear,
                make sure that ``my_function(0)`` is finite.

            .. hint::
                In order to mimic a singular user defined function,
                use the parameter `regulator_power`.
                For example, instead of ``numerator = gamma(eps)`` you
                could enter ``numerator = eps_times_gamma(eps)`` in
                conjunction with ``regulator_power = -1``

            .. warning::
                The `numerator` is very flexible in its input. However,
                that flexibility comes for the price of less error
                safety. We have no way of checking for all possible
                mistakes in the input. If your numerator is more
                advanced than in the examples above, you should proceed
                with great caution.

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

            .. note::
                If you change this symbol, you have to adapt
                the `dimensionality` accordingly.

        :param regulator_power:
            integer;
            The regulator to this power will be multiplied by
            the numerator.

            .. seealso::
                parameter `numerator`

        :param dimensionality:
            string or sympy expression, optional;
            The dimensionality; typically :math:`4-2\epsilon`,
            which is the default value.

        :param metric_tensor:
            string or sympy symbol, optional;
            The symbol to be used for the (Minkowski) metric
            tensor :math:`g^{\mu\nu}`.

    '''
        # TODO: carefully reread and check this documentation
        # TODO: test case with linear propagators
        # TODO: test that including all the factors described above
        #       is the complete integrand

    def __init__(self, propagators, loop_momenta, external_momenta=[], Lorentz_indices=[], \
                 numerator=1, replacement_rules=[], Feynman_parameters='x', regulator='eps', \
                 regulator_power=0, dimensionality='4-2*eps', metric_tensor='g'):

        def sympify_symbols(iterable, error_message, allow_number=False):
            '''
            `sympify` each item in `iterable` and assert
            that it is a `symbol`.

            '''
            symbols = []
            for expression in iterable:
                expression = sp.sympify(expression)
                assert expression.is_Symbol or expression.is_Number if allow_number else expression.is_Symbol, error_message
                symbols.append(expression)
            return symbols

        def assert_at_most_quadractic(expression, variables, error_message):
            '''
            Assert that `expression` is a polynomial of
            degree less or equal 2 in the `variables`.

            '''
            poly = Polynomial.from_expression(expression, variables)
            assert (poly.expolist.sum(axis=1) <= 2).all(), error_message

        # check and store `regulator_power`
        regulator_power_as_int = int(regulator_power)
        assert regulator_power_as_int == regulator_power, '`regulator_power` must be integral.'
        self.regulator_power = regulator_power_as_int

        # sympify and store `loop_momenta`
        self.loop_momenta = sympify_symbols(loop_momenta, 'Each of the `loop_momenta` must be a symbol.')
        self.L = len(self.loop_momenta)

        # sympify and store `external_momenta`
        self.external_momenta = sympify_symbols(external_momenta, 'Each of the `external_momenta` must be a symbol.')

        all_momenta = self.external_momenta + self.loop_momenta

        # sympify and store `Lorentz_indices`
        self.Lorentz_indices = sympify_symbols(Lorentz_indices, 'Each of the `Lorentz_indices` must be a symbol or a number.', allow_number=True)

        # sympify and store `metric_tensor`
        self.metric_tensor = sympify_symbols([metric_tensor], '`metric_tensor` must be a symbol.')[0]

        # sympify and store `regulator`
        self.regulator = sympify_symbols([regulator], '`regulator` must be a symbol.')[0]

        # sympify and store `dimensionality`
        self.dimensionality = sp.sympify(dimensionality)

        # sympify and store `propagators`
        self.propagators = sp.sympify(list(propagators))
        for propagator in self.propagators:
            assert_at_most_quadractic(propagator, all_momenta, 'Each of the `propagators` must be polynomial and at most quadratic in the momenta.')
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

        if not isinstance(replacement_rules, list):
            replacement_rules = list(replacement_rules)
        if replacement_rules:
            self.replacement_rules = np.array(replacement_rules)
            assert len(self.replacement_rules.shape) == 2, "The `replacement_rules` should be a list of tuples"
            assert self.replacement_rules.shape[1] == 2 , "The `replacement_rules` should be a list of tuples"
            for rule in self.replacement_rules:
                for expression in rule:
                    assert_at_most_quadractic(expression, all_momenta, 'Each of the `replacement_rules` must be polynomial and at most quadratic in the momenta.')
        else:
            self.replacement_rules = []

        self.numerator_input = sp.sympify(numerator).expand()
        if self.numerator_input.is_Add:
            self.numerator_input_terms = list(self.numerator_input.args)
        else:
            self.numerator_input_terms = [self.numerator_input]


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
                # all coeffs of the polynomials in M are likely to be integer
                #    --> try conversion to `int` since it is faster than calculating with sympy expressions
                try:
                    tmp.coeffs = tmp.coeffs.astype(int)
                # ignore if not all coefficients can be converted to `int` --> keep the sympy expressions
                except:
                    pass
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
        def append_double_indices(expr, lst, momenta, indices):
            '''
            Append the double indices of an expression of the form
            ``<momentum>(<index>)`` to `lst`.
            Return the remaining factors.

            '''
            for index in indices:
                index_count = 0
                for i,momentum in enumerate(momenta):
                    for j in range(2): # should have the same index at most twice
                        if expr.subs(momentum(index), 0) == 0: # expression has a factor `momentum(index)`
                            expr /= momentum(index)
                            lst.append((i,index))
                            index_count += 1
                        else:
                            break
                    # should not have factors of `momentum(index)` left after the loop above
                    assert str(momentum(index)) not in str(expr), 'Could not parse `numerator`'
                assert index_count <= 2, 'Could not parse `numerator`'
            return expr

        numerator_loop_tensors = []
        numerator_external_tensors = []
        numerator_symbols = []
        for term in self.numerator_input_terms:
            current_loop_tensor = []
            current_external_tensor = []

            # search for ``momentum(index)``
            term = append_double_indices(term, current_loop_tensor, self.loop_momenta, self.Lorentz_indices)
            term = append_double_indices(term, current_external_tensor, self.external_momenta, self.Lorentz_indices)
            numerator_loop_tensors.append(current_loop_tensor)
            numerator_external_tensors.append(current_external_tensor)
            numerator_symbols.append(term)
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
    def numerator_ranks(self):
        return [len(indices) for indices in self.numerator_loop_tensors]

    @cached_property
    def highest_rank(self):
        return max(self.numerator_ranks)

    @cached_property
    def numerator(self):
        '''
        Generate the numerator in index notation according to
        section 2 in [GKR+11]_.

        '''
        # The implementation below conceptually follows section 2 of arXiv:1010.1667.

        numerator = 0

        aM = self.aM
        Q = self.Q
        g = self.metric_tensor
        D = self.dimensionality
        L = self.L
        Feynman_parameters_F_U = self.Feynman_parameters + ['F', 'U']
        U = Polynomial.from_expression('U', Feynman_parameters_F_U)
        F = Polynomial.from_expression('F', Feynman_parameters_F_U)
        replacement_rules = self.replacement_rules_with_Lorentz_indices
        highest_rank = self.highest_rank
        N_nu = len(self.propagators)

        # Every term factor in the sum of equation (2.5) in arXiv:1010.1667v1 comes with
        # the scalar factor `1/(-2)**(r/2)*Gamma(N_nu - D*L/2 - r/2)*F**(r/2)`.
        # In order to keep the `numerator` free of poles in the regulator, we divide it
        # by the Gamma function with the smallest argument `N_nu - D*L/2 - highest_rank//2`,
        # where `//` means integer division, and use `x*Gamma(x) = Gamma(x+1)`.
        one_over_minus_two = sp.sympify('1/(-2)')
        def scalar_factor(r):
            # `r` must be even
            assert r % 2 == 0
            factors_without_gamma = one_over_minus_two**(r//2) * F**(r//2)
            # gamma_factor = 'Gamma(N_nu - D*L/2 - r/2)'
            def reduce_gamma(r_over_two):
                '''
                   Transform Gamma(N_nu - D*L/2 - r_over_two) -->
                   (N_nu - D*L/2 - (r_over_two+1))*(N_nu - D*L/2 - (r_over_two+2))*...*Gamma(N_nu - D*L/2 - highest_rank//2)'
                   and divide by Gamma(N_nu - D*L/2 - highest_rank//2).
                '''
                if r_over_two == highest_rank//2:
                    return 1
                return (N_nu - D*L/2 - (r_over_two+1)) * reduce_gamma(r_over_two+1)
            return factors_without_gamma * reduce_gamma(r//2)

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

            # Since we allow sums of numerators with different rank, we must multiply the correct
            # power of U. We globally multiply by ``U**(N_nu - dim / 2 * (L + 1) - highest_rank)``
            # Therefore we must multiply each term by ``U**(highest_rank - this_term_rank)``
            this_U_factor = U**( highest_rank - len(loop_tensor) )

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

                    this_numerator_summand = scalar_factor(r) * this_U_factor * remainder

                    # --------------------------- multiply the tensor "P" --------------------------
                    for external_momentum_index, Lorentz_index in P_indices:
                        this_tensor_P_factor = aM[external_momentum_index].dot(Q)
                        # There should be eactly one external momentum in each term (summand) of `this_tensor_P_factor` --> attach the `Lorentz_index` to it
                        for i, coeff in enumerate(this_tensor_P_factor.coeffs):
                            this_tensor_P_factor.coeffs[i] = coeff.subs((p, p(Lorentz_index)) for p in self.external_momenta)

                        # must append ``F`` and ``U`` to the parameters of ``this_tensor_P_factor``
                        this_tensor_P_factor.expolist = np.hstack([this_tensor_P_factor.expolist, np.zeros((len(this_tensor_P_factor.expolist), 2), dtype=int)])
                        this_tensor_P_factor.number_of_variables += 2

                        this_numerator_summand *= this_tensor_P_factor
                    # ------------------------------------------------------------------------------

                    # ---- multiply the caligraphic "A" in equation (2.5) of arXiv:1010.1667v1. ----
                    # `tensors_A2_indices` is a list of pairs of double indices describing one of the factors of equation (2.15) in arXiv:1010.1667v1.
                    # The following loop inserts the definition (2.18, arXiv:1010.1667v1) for all of these pairs.
                    # The metric tensors in "A" can contract momenta of the `external_tensor` --> generate a list of
                    # contractions and multiply only the uncontracted metric tensors
                    contractions_left = []
                    contractions_right = []
                    contractions_unmatched = []
                    for tensor_A2_indices in tensors_A2_indices:
                        aM_factor = aM[tensor_A2_indices[0][0],tensor_A2_indices[1][0]] # * g(tensor_A2_indices[0][1], tensor_A2_indices[1][1])

                        # must append the variables ``F`` and ``U`` to ``aM_factor``
                        if isinstance(aM_factor, Polynomial):
                            aM_factor = aM_factor.copy()
                            aM_factor.expolist = np.hstack([aM_factor.expolist, np.zeros((len(aM_factor.expolist), 2), dtype=int)])
                            aM_factor.number_of_variables += 2

                        this_numerator_summand *= aM_factor
                        contractions_left.append(tensor_A2_indices[0][1])
                        contractions_right.append(tensor_A2_indices[1][1])
                        contractions_unmatched.append(True)
                    # ------------------------------------------------------------------------------

                    # ----------------------- multiply the `external_tensor` -----------------------
                    for external_momentum_index, Lorentz_index_external in external_tensor:
                        # replace metric tensors where possible
                        external_momentum_unmatched = True
                        for i, (Lorentz_index_1, Lorentz_index_2, metric_unmatched) in enumerate(zip(contractions_left, contractions_right, contractions_unmatched)):
                            if metric_unmatched:
                                if Lorentz_index_1 == Lorentz_index_external:
                                    contractions_unmatched[i] = external_momentum_unmatched = False
                                    this_numerator_summand *= self.external_momenta[external_momentum_index](Lorentz_index_2)
                                    break
                                if Lorentz_index_2 == Lorentz_index_external:
                                    contractions_unmatched[i] = external_momentum_unmatched = False
                                    this_numerator_summand *= self.external_momenta[external_momentum_index](Lorentz_index_1)
                                    break
                        if external_momentum_unmatched:
                            this_numerator_summand *= self.external_momenta[external_momentum_index](Lorentz_index_external)
                    # ------------------------------------------------------------------------------

                    # -------------------------- multiply "A" (continued) --------------------------
                    # multiply all unmatched metric tensors
                    for Lorentz_index_1, Lorentz_index_2, is_unmatched in zip(contractions_left, contractions_right, contractions_unmatched):
                        if is_unmatched:
                            if Lorentz_index_1 == Lorentz_index_2:
                                # g(mu, mu) = dimensionality
                                this_numerator_summand *= D
                            else:
                                this_numerator_summand *= g(Lorentz_index_1, Lorentz_index_2)
                    # ------------------------------------------------------------------------------

                    # apply the replacement rules
                    for i, coeff in enumerate(this_numerator_summand.coeffs):
                        this_numerator_summand.coeffs[i] = coeff.expand().subs(replacement_rules)

                    numerator += this_numerator_summand

        # The global factor of `(-1)**N_nu` in equation (2.15) of arXiv:1010.1667v1 must appear somewhere.
        # The numerator seems a good choice since this has a nonsingular but nontrivial expansion in the
        # regulator for noninteger propagator powers.
        return numerator * ( (-1)**N_nu )

    @cached_property
    def replacement_rules_with_Lorentz_indices(self):
        replacement_rules = []
        for rule in self.replacement_rules:
            pattern = sp.sympify(rule[0])
            replacement = sp.sympify(rule[1])
            for mu in self.Lorentz_indices:
                pattern_with_index = pattern
                replacement_with_index = replacement
                for p in self.external_momenta:
                    pattern_with_index = pattern_with_index.subs(p, p(mu))
                    replacement_with_index = replacement_with_index.subs(p, p(mu))
                replacement_rules.append( (pattern_with_index, replacement_with_index) )

        return replacement_rules


class LoopIntegral_from_graph(LoopIntegral):
    '''
    Construct the Feynman parametrization of a
    loop integral from the graph using the cut construction method.
    '''
    # TODO: implement replacement rules in cut construct

    def __init__(self, loops, lines, external_momenta):
    
        self.L=loops

        self.intlines=[]
        for line in lines:
            temp=[sp.sympify(line[0]),[line[1],line[2]]]
            self.intlines.append(temp)

        self.extlines=[]
        for i in range(len(external_momenta)):
            temp=[sp.sympify(external_momenta[i][0]),[-i-1,external_momenta[i][1]]]
            self.extlines.append(temp)

    @cached_property
    def vertdict(self): 
        # creates a list of internal and external vertices and indexes them
        lines = self.intlines + self.extlines
        vertices=set([])
        for line in lines:
            vertices = vertices | set(line[1])
        vertices=list(vertices)
        dict = {}
        for i in range(len(vertices)):
            dict[vertices[i]] = i
        return dict
        
    @cached_property
    def vertmatrix(self):
        # create transition matrix representation of underlying graph
        numvert = len(self.vertdict)
        M = np.empty((numvert,numvert), dtype=object)
        # each vertex is trivially connected to itself
        for i in range(numvert):
            for j in range(i,numvert):
                if i == j:
                    M[i,j]=1
                else:
                    M[i,j]=M[j,i]=0
        for i in range(len(self.intlines)):
            # get indices for start and end of each propagator
            start = self.vertdict[self.intlines[i][1][0]]
            end = self.vertdict[self.intlines[i][1][1]]
            # add entry to the matrix representation
            temp = sp.sympify('pr'+str(i))
            M[start,end] += temp
            M[end,start] += temp
        return M
        
    @cached_property
    def U(self):
        raise NotImplementedError()
        # TODO: implement cut construction

    @cached_property
    def F(self):
        raise NotImplementedError()
        # TODO: implement cut construction
        
