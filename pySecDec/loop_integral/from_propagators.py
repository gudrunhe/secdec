"""Routines to Feynman parametrize a loop integral given its algebraic propagators."""

from .common import LoopIntegral
from ..algebra import Polynomial
from ..misc import det, adjugate, powerset, missing, all_pairs, \
     cached_property, sympify_symbols, assert_degree_at_most_max_degree, sympify_expression, rec_subs
import sympy as sp
import numpy as np

# sympy symbols are no longer callable starting from version 1.3
_to_function = lambda x: sp.Function(str(x))

class LoopIntegralFromPropagators(LoopIntegral):
    __doc__ = r'''
    Construct the Feynman parametrization of a
    loop integral from the algebraic momentum
    representation.

    .. seealso::
        [Hei08]_, [GKR+11]_

    Example:

    >>> from pySecDec.loop_integral import *
    >>> propagators = ['k**2', '(k - p)**2']
    >>> loop_momenta = ['k']
    >>> li = LoopIntegralFromPropagators(propagators, loop_momenta)
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
        The propagators, e.g. ``['k1**2', '(k1-k2)**2 - m1**2']``.

    :param loop_momenta:
        iterable of strings or sympy expressions;
        The loop momenta, e.g. ``['k1','k2']``.

    :param external_momenta:
        iterable of strings or sympy expressions,
        optional;
        The external momenta, e.g. ``['p1','p2']``.
        Specifying the `external_momenta` is only
        required when a `numerator` is to be
        constructed.

        .. seealso::
            parameter `numerator`

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
        ``k1(mu)*k2(mu)``. The numerator must be a sum of products
        of exclusively:
         * numbers
         * scalar products (e.g. ``p1(mu)*k1(mu)*p1(nu)*k2(nu)``)
         * symbols (e.g. ``s``, ``eps``)

        Examples:
            * ``p1(mu)*k1(mu)*p1(nu)*k2(nu) + 4*s*eps*k1(mu)*k1(mu)``
            * ``p1(mu)*(k1(mu) + k2(mu))*p1(nu)*k2(nu)``
            * ``p1(mu)*k1(mu)``

        .. note::
            In order to use the resulting `LoopIntegral` as an argument
            to the function :func:`pySecDec.loop_integral.loop_package`,
            the resulting Feynman parametrized ``self.numerator`` must
            be expressible as a :class:`pySecDec.algebra.Polynomial`
            such that all coefficients are purely numeric. In
            addition, all scalar products of the external momenta
            must be expressed in terms of Mandelstam variables using
            the `replacement_rules`.

        .. warning::
            **All** Lorentz indices (including the contracted ones
            and also including the numbers that have been used)
            must be explicitly defined using the parameter
            `Lorentz_indices`.

        .. hint::
            It is possible to use numbers as indices, for example
            ``p1(mu)*p2(mu)*k1(nu)*k2(nu) = p1(1)*p2(1)*k1(2)*k2(2)``.

        .. hint::
            The numerator may have uncontracted indices, e.g.
            ``k1(mu)*k2(nu)``. If indices are left open, however, the
            LoopIntegral cannot be used with the package generator
            :func:`pySecDec.loop_integral.loop_package`.

    :param metric_tensor:
        string or sympy symbol, optional;
        The symbol to be used for the (Minkowski) metric
        tensor :math:`g^{\mu\nu}`.

    ''' + LoopIntegral.common_properties_doc

    def __init__(self, propagators, loop_momenta, external_momenta=[], Lorentz_indices=[],
                 numerator=1, metric_tensor='g', replacement_rules=[], Feynman_parameters='x',
                 regulators=None, regulator=None,
                 dimensionality='4-2*eps', powerlist=[]):

        # sympify and store `loop_momenta`
        self.loop_momenta = sympify_symbols(loop_momenta, 'Each of the `loop_momenta` must be a symbol.')
        self.L = len(self.loop_momenta)

        # sympify and store `external_momenta`
        self.external_momenta = sympify_symbols(external_momenta, 'Each of the `external_momenta` must be a symbol.')

        self.all_momenta = self.external_momenta + self.loop_momenta

        # sympify and store `Lorentz_indices`
        self.Lorentz_indices = sympify_symbols(Lorentz_indices, 'Each of the `Lorentz_indices` must be a symbol or a number.', allow_number=True)

        # sympify and store `metric_tensor`
        self.metric_tensor = sympify_symbols([metric_tensor], '`metric_tensor` must be a symbol.')[0]

        # sympify and store `propagators`
        self.propagators = sympify_expression(list(propagators))
        for propagator in self.propagators:
            assert_degree_at_most_max_degree(propagator, self.all_momenta, 2, 'Each of the `propagators` must be polynomial and at most quadratic in the momenta.')
        self.P = len(self.propagators)

        # sympify and store `numerator`
        self.numerator_input = sympify_expression(numerator).expand()
        if self.numerator_input.is_Add:
            self.numerator_input_terms = list(self.numerator_input.args)
        else:
            self.numerator_input_terms = [self.numerator_input]

        # store properties shared between derived classes
        regulators = \
                regulators if regulators is not None else \
                [regulator] if regulator is not None else \
                ["eps"]
        self.set_common_properties(replacement_rules, Feynman_parameters, regulators,
                                   dimensionality, powerlist)

        # remove `propagators` and `Feynman_parameters` that are set to zero by the `powerlist`
        for i in range(self.P-1,-1,-1): # traverse backwards to stay consistent with the indexing
            if self.powerlist[i] == 0:
                self.P -= 1
                self.powerlist = self.powerlist[:i] + self.powerlist[i+1:]
                self.propagators = self.propagators[:i] + self.propagators[i+1:]
                self.derivativelist = self.derivativelist[:i] + self.derivativelist[i+1:]
                self.Feynman_parameters = self.Feynman_parameters[:i] + self.Feynman_parameters[i+1:]

        # If the degree of an inverse propagator is smaller than 2,
        # the derivative of U with respect to the corresponding Feynman parameter will vanish,
        # which will change the exponent of U.
        for i in range(self.P):
            k = self.derivativelist[i]
            if k !=0 :
                poly = Polynomial.from_expression(self.propagators[i], self.loop_momenta)
                if max(poly.expolist.sum(axis=1)) < 2:
                    self.U_derivatives -= k

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
    preliminary_U = detM

    @cached_property
    def preliminary_F(self):
        # equation (8) of arXiv:0803.4177: F = det(M)*(Q.transpose*inverse(M)*Q-J) = (Q.transpose*adjugate(M)*Q-U*J)
        F = 0
        for i in range(self.L):
            for j in range(self.L):
                F += self.Q[i]*self.aM[i,j]*self.Q[j]
        F -= self.preliminary_U * self.J
        for i,coeff in enumerate(F.coeffs):
            if isinstance(coeff, sp.Expr):
                F.coeffs[i] = rec_subs(coeff,self.replacement_rules)
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
        numerator_loop_tensors = []
        numerator_external_tensors = []
        numerator_symbols = []
        for term in self.numerator_input_terms:
            original_term = term
            current_loop_tensor = []
            current_external_tensor = []

            for index in self.Lorentz_indices:
                index_count = 0
                for lst,momenta in zip([current_loop_tensor,current_external_tensor],[self.loop_momenta, self.external_momenta]):
                    for i,momentum in enumerate(momenta):

                        # search for ``momentum(index)``
                        while term.subs(_to_function(momentum)(index), 0) == 0: # expression has a factor `momentum(index)`
                            term /= _to_function(momentum)(index)
                            lst.append((i,index))
                            index_count += 1

                            # should not have found the `index` more than twice
                            assert index_count <= 2, \
                                'A term (%s) in the `numerator` has one of the `Lorentz_indices` (%s) occurring more than twice.' % (original_term, index)

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
    def preliminary_numerator(self):
        '''
        Generate the numerator in index notation according to
        section 2 in [GKR+11]_.

        '''
        # The implementation below conceptually follows section 2 of arXiv:1010.1667.

        numerator = 0

        aM = self.aM
        Q = self.Q
        g = _to_function(self.metric_tensor)
        D = self.dimensionality
        L = self.L
        Feynman_parameters_U_F = self.Feynman_parameters + sympify_expression(['U', 'F'])
        U = Polynomial.from_expression('U', Feynman_parameters_U_F)
        F = Polynomial.from_expression('F', Feynman_parameters_U_F)
        replacement_rules = self.replacement_rules_with_Lorentz_indices
        highest_rank = self.highest_rank
        N_nu = self.N_nu

        if self.numerator_input == 1:
            return Polynomial(np.zeros([1,len(Feynman_parameters_U_F)], dtype=int), np.array([1]), Feynman_parameters_U_F, copy=False)

        # Every term factor in the sum of equation (2.5) in arXiv:1010.1667v1 comes with
        # the scalar factor `1/(-2)**(r/2)*Gamma(N_nu - D*L/2 - r/2)*F**(r/2)`.
        # In order to keep the `numerator` free of poles in the regulator, we divide it
        # by the Gamma function with the smallest argument `N_nu - D*L/2 - highest_rank//2`,
        # where `//` means integer division, and use `x*Gamma(x) = Gamma(x+1)`.
        one_over_minus_two = sympify_expression('1/(-2)')
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

                    # match ``g(mu, nu) * g(rho, sigma)`` when two indices in different ``g``s are equal
                    matched = True
                    while matched:
                        matched = False

                        for i, (Lorentz_index_mu, Lorentz_index_nu, left_metric_unmatched) in enumerate(zip(contractions_left, contractions_right, contractions_unmatched)):
                            if not left_metric_unmatched:
                                continue

                            for j, (Lorentz_index_rho, Lorentz_index_sigma, right_metric_unmatched) in enumerate(zip(contractions_left, contractions_right, contractions_unmatched)):

                                # we must have ``i > j`` because python seems to not take e.g. ``contractions_left[j] = ...`` into account in the iterator
                                if i <= j or not right_metric_unmatched:
                                    continue

                                # ``mu == rho`` -- > ``g(mu, nu) * g(mu, sigma) = g(nu, sigma)``
                                if Lorentz_index_mu == Lorentz_index_rho:
                                    # re-execute while as long as there are simplifications
                                    matched = True
                                    # remove ``g(mu, nu)``
                                    contractions_unmatched[i] = False
                                    # "Lorentz_index_rho --> Lorentz_index_nu"; i.e. replace ``g(rho, sigma) --> g(nu, sigma)``
                                    contractions_left[j] = Lorentz_index_nu
                                    break

                                # ``mu == sigma`` --> ``g(mu, nu) * g(rho, mu) = g(rho, nu)``
                                elif Lorentz_index_mu == Lorentz_index_sigma:
                                    # re-execute while as long as there are simplifications
                                    matched = True
                                    # remove ``g(mu, nu)``
                                    contractions_unmatched[i] = False
                                    # "Lorentz_index_sigma --> Lorentz_index_nu"; i.e. replace ``g(rho, sigma) --> g(rho, nu)``
                                    contractions_right[j] = Lorentz_index_nu
                                    break

                                # ``nu == rho`` --> ``g(mu, nu) * g(nu, sigma) = g(mu, sigma)``
                                elif Lorentz_index_nu == Lorentz_index_rho:
                                    # re-execute while as long as there are simplifications
                                    matched = True
                                    # remove ``g(mu, nu)``
                                    contractions_unmatched[i] = False
                                    # "Lorentz_index_rho --> Lorentz_index_mu"; i.e. replace ``g(rho, sigma) --> g(mu, sigma)``
                                    contractions_left[j] = Lorentz_index_mu
                                    break

                                # ``nu == sigma`` --> ``g(mu, nu) * g(rho, nu) = g(rho, mu)``
                                elif Lorentz_index_nu == Lorentz_index_sigma:
                                    # re-execute while as long as there are simplifications
                                    matched = True
                                    # remove ``g(mu, nu)``
                                    contractions_unmatched[i] = False
                                    # "Lorentz_index_sigma --> Lorentz_index_mu"; i.e. replace ``g(rho, sigma) --> g(rho, mu)``
                                    contractions_right[j] = Lorentz_index_mu
                                    break

                    # --------------------------- multiply the tensor "P" --------------------------

                    # There should be exactly one external momentum in each term (summand) of `this_tensor_P_factor`

                    # --> contract with metric tensors
                    # we need a copy of the ``P_indices``
                    contracted_P_indices = list(P_indices)
                    for k, (external_momentum_index, Lorentz_index_P) in enumerate(P_indices):
                        for i, (Lorentz_index_1, Lorentz_index_2, metric_unmatched) in enumerate(zip(contractions_left, contractions_right, contractions_unmatched)):
                            if metric_unmatched:
                                if Lorentz_index_1 == Lorentz_index_P:
                                    contractions_unmatched[i] = False
                                    contracted_P_indices[k] = (external_momentum_index, Lorentz_index_2)
                                    break
                                if Lorentz_index_2 == Lorentz_index_P:
                                    contractions_unmatched[i] = False
                                    contracted_P_indices[k] = (external_momentum_index, Lorentz_index_1)
                                    break

                    # --> attach the `Lorentz_index_P` to it
                    for k, (external_momentum_index, Lorentz_index_P) in enumerate(contracted_P_indices):
                        this_tensor_P_factor = aM[external_momentum_index].dot(Q)
                        for j, coeff in enumerate(this_tensor_P_factor.coeffs):
                            if isinstance(coeff, sp.Expr):
                                this_tensor_P_factor.coeffs[j] = coeff.subs((p, _to_function(p)(Lorentz_index_P)) for p in self.external_momenta)

                        # must append ``F`` and ``U`` to the parameters of ``this_tensor_P_factor``
                        this_tensor_P_factor.expolist = np.hstack([this_tensor_P_factor.expolist, np.zeros((len(this_tensor_P_factor.expolist), 2), dtype=int)])
                        this_tensor_P_factor.number_of_variables += 2

                        this_numerator_summand *= this_tensor_P_factor
                    # ------------------------------------------------------------------------------

                    # ----------------------- multiply the `external_tensor` -----------------------

                    # replace metric tensors where possible

                    # match ``g(mu, nu)`` with external_momenta (e.g. ``p(mu)``)
                    for external_momentum_index, Lorentz_index_external in external_tensor:
                        external_momentum_unmatched = True
                        for i, (Lorentz_index_1, Lorentz_index_2, metric_unmatched) in enumerate(zip(contractions_left, contractions_right, contractions_unmatched)):
                            if metric_unmatched:
                                if Lorentz_index_1 == Lorentz_index_external:
                                    contractions_unmatched[i] = external_momentum_unmatched = False
                                    this_numerator_summand *= _to_function(self.external_momenta[external_momentum_index])(Lorentz_index_2)
                                    break
                                if Lorentz_index_2 == Lorentz_index_external:
                                    contractions_unmatched[i] = external_momentum_unmatched = False
                                    this_numerator_summand *= _to_function(self.external_momenta[external_momentum_index])(Lorentz_index_1)
                                    break
                        if external_momentum_unmatched:
                            this_numerator_summand *= _to_function(self.external_momenta[external_momentum_index])(Lorentz_index_external)

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
                        this_numerator_summand.coeffs[i] = rec_subs(sympify_expression(coeff), replacement_rules)

                    numerator += this_numerator_summand

        return numerator

    @cached_property
    def replacement_rules_with_Lorentz_indices(self):
        replacement_rules = []
        for rule in self.replacement_rules:
            pattern = sympify_expression(rule[0])
            replacement = sympify_expression(rule[1])
            for mu in self.Lorentz_indices:
                pattern_with_index = pattern
                replacement_with_index = replacement
                for p in self.external_momenta:
                    pattern_with_index = pattern_with_index.subs(p, _to_function(p)(mu))
                    replacement_with_index = replacement_with_index.subs(p, _to_function(p)(mu))
                replacement_rules.append( (pattern_with_index, replacement_with_index) )

        return replacement_rules
