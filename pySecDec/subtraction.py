"""
Subtraction
-----------

Routines to isolate the divergencies in an :math:`\epsilon`
expansion.

"""

from .algebra import Polynomial, ExponentiatedPolynomial, Sum, Product, Pow
from .misc import sympify_expression
from itertools import repeat
import numpy as np

_sympy_one = sympify_expression(1)

def pole_structure(monomial_product, *indices):
    '''
    Return a list of the unregulated exponents
    of the parameters specified by `indices`
    in `monomial_product`.

    :param monomial_product:
        :class:`pySecDec.algebra.ExponentiatedPolynomial`
        with ``exponent`` being a :class:`.Polynomial`;
        The monomials of the subtraction to extract the
        pole structure from.

    :param indices:
        arbitrarily many integers;
        The index/indices of the parameter(s) to partially
        investigate.

    '''
    def pole_structure_one_index(monomial_product, index):
        exponent_constant_term = 0
        for monomial_factor in monomial_product.factors:
            if monomial_factor.exponent.has_constant_term():
                this_factor_exponent_constant_term = np.sum(monomial_factor.exponent.coeffs[np.where((monomial_factor.exponent.expolist == 0).all(axis=1))])
                exponent_constant_term += this_factor_exponent_constant_term * monomial_factor.expolist[0,index]
            # else: e.g. (z1**2)**-eps0-eps1 --> no contribution to ``a_j``
                # exponent_constant_term += 0
        return exponent_constant_term

    return [pole_structure_one_index(monomial_product, index) for index in indices]

# -------------------------------------------- original subtraction --------------------------------------------

def _integrate_pole_part_single_index(polyprod, index):
    monomial_product = polyprod.factors[0] # e.g. (z1**2)**1-eps0-eps1  *  (z1*z2)**3-4*eps1
    monomial_product_FeynmanJ_set_to_one = monomial_product.replace(index,1)
    regulator_poles = polyprod.factors[1]
    cal_I = Product(*polyprod.factors[2:], copy=False)
    polysymbols = monomial_product.symbols

    # arXiv:0803.4177v2: exponent_constant_term = a_j
    exponent_constant_term = 0
    full_exponent = 0
    for monomial_factor in monomial_product.factors:
        if monomial_factor.exponent.has_constant_term():
            this_factor_exponent_constant_term = np.sum(monomial_factor.exponent.coeffs[np.where((monomial_factor.exponent.expolist == 0).all(axis=1))])
            exponent_constant_term += this_factor_exponent_constant_term * monomial_factor.expolist[0,index]
        # else: e.g. (z1**2)**-eps0-eps1 --> no contribution to ``a_j``
            # exponent_constant_term += 0
        full_exponent += monomial_factor.exponent * monomial_factor.expolist[0,index]

    if exponent_constant_term < 0:
        regulator_part = full_exponent - exponent_constant_term
        # the next line conceptually implements "assert regulator_part != 0"
        if (regulator_part.coeffs == 0).all():
            raise ValueError('"1/0" detected. An additional regulator is needed to calculate this integral.')

    # finite variance of Monte Carlo integral estimator only if ``exponent_constant_term > -0.5``
    # to be absolutely safe here, we eliminate integrable singularities; i.e. ``exponent_constant_term >= 0``
    # TODO: Change to a dynamic `power_goal` just like in :func:`.integrate_by_parts`
    if exponent_constant_term >= 0:
        # no subtraction needed, the input `polyprod` is numerically integrable
        return [polyprod]

    def make_FeynmanIndex_to_power(power, factorial_of_power, symbols):
        '''
        Return a :class:`Polynomial` representing
        "Feynman_index**power/power!"

        '''
        expolist = np.zeros((1, monomial_product.number_of_variables), dtype=int)
        expolist[:,index] = power
        return Polynomial(expolist, np.array([_sympy_one/factorial_of_power]), symbols, copy=False)

    output_summands = []
    minus_cal_I_expansion_summands = []
    # exponent_constant_term <= -1 < 0      =>     abs(exponent_constant_term) = -exponent_constant_term
    # use symbol names as in arXiv:0803.4177v2

    # construction of the pole part
    for p in range(int(-exponent_constant_term)):
        # No dependency on Feynman parameter with index `index` --> has been integrated out analytically
        #   --> Replace that Feynman parameter by `1` in the monomial factors.
        current_factors = [monomial_product_FeynmanJ_set_to_one]

        # renew `p_factorial` on the fly
        p_factorial = 1 if p == 0 else p_factorial * p

        # renew derivative of cal_I on the fly
        derivative_cal_I = cal_I if p == 0 else derivative_cal_I.derive(index)
        derivative_cal_I_Feynmanj_set_to_zero = derivative_cal_I.replace(index, 0)

        # arXiv0803.4177v2: 1/( (a_j + p + 1 - b_j * eps) * factorial(p) )
        new_potential_pole_denominator = (full_exponent + (p + 1)) * (p_factorial)
        # put this factor into the pole part only if a_j + p + 1 is zero
        if exponent_constant_term + p + 1 == 0:
            current_regulator_poles = Pow(regulator_poles.base * new_potential_pole_denominator, regulator_poles.exponent, copy=False)
            current_factors.append(current_regulator_poles)
            current_factors.append(derivative_cal_I_Feynmanj_set_to_zero.copy())
        # otherwise it does not lead to additional regulator poles and can become part of <derivative_cal_I>
        else:
            # poles in current term: none --> poles are just the old ones
            current_factors.append(regulator_poles)

            # put `new_potential_pole_denominator**-1` (which is not a pole in this case) in the last factor
            last_factor = Product(derivative_cal_I_Feynmanj_set_to_zero.copy(), Pow(new_potential_pole_denominator, regulator_poles.exponent.copy(), copy=False), copy=False)
            current_factors.append(last_factor.simplify())

        output_summands.append(Product(*current_factors, copy=False))
        minus_cal_I_expansion_summands.append(Product(-make_FeynmanIndex_to_power(p, p_factorial, polysymbols), derivative_cal_I_Feynmanj_set_to_zero, copy=False ))

    integrable_part = Product(   monomial_product, regulator_poles, Sum(cal_I, *minus_cal_I_expansion_summands, copy=False).simplify() , copy=False   )

    output_summands.append(integrable_part)

    return output_summands

def integrate_pole_part(polyprod, *indices):
    r'''
    Transform an integral of the form

    .. math::
        \int_0^1
        {
            dt_j t_j^{(a - b \epsilon_1 - c \epsilon_2 + ...)}
            \mathcal{I} (t_j,\{t_{i \neq j}\}, \epsilon_1, \epsilon_2, ...)
        }

    into the form

    .. math::
        \sum_{p=0}^{|a|-1}
        {
            \frac{1}{a + p + 1 - b \epsilon_1 - c \epsilon_2 - ...}
            \frac{\mathcal{I}^{(p)} (0,\{t_{i \neq j}\}, \epsilon_1, \epsilon_2, ...)}{p!} +
            \int_0^1
            {
                dt_j t_j^{(a - b \epsilon_1 - c \epsilon_2 + ...)}
                R(t_j,\{t_{i \neq j}\}, \epsilon_1, \epsilon_2, ...)
            }
        }

    , where :math:`\mathcal{I}^{(p)}` denotes the p-th derivative
    of :math:`\mathcal{I}` with respect to :math:`t_j`. The equations
    above are to be understood schematically.

    .. seealso::
        This function implements the transformation from
        equation (19) to (21) as described in arXiv:0803.4177v2
        [Hei08]_.


    :param polyprod:
        :class:`.algebra.Product` of the
        form ``<product of <monomial>**(a_j + ...)> *
        <regulator poles of cal_I> * <cal_I>``;
        The input product as decribed above.
        The <product of <monomial>**(a_j + ...)> should be
        a :class:`pySecDec.algebra.Product` of
        <monomial>**(a_j + ...).
        as described below.
        The <monomial>**(a_j + ...) should be an
        :class:`pySecDec.algebra.ExponentiatedPolynomial`
        with ``exponent`` being a :class:`.Polynomial` of the
        regulators :math:`\epsilon_1, \epsilon_2, ...`. Although
        no dependence on the Feynman parameters is expected
        in the ``exponent``, the polynomial variables should
        be the Feynman parameters and the regulators.
        The constant term of the exponent should be numerical.
        The polynomial variables of ``monomial`` and the other
        factors (interpreted as :math:`\mathcal{I}`) are interpreted
        as the Feynman parameters and the epsilon regulators.
        Make sure that the last factor (``<cal_I>``) is defined
        and finite for :math:`\epsilon = 0`. All poles for
        :math:`\epsilon \rightarrow 0` should be made explicit
        by putting them into ``<regulator poles of cal_I>``
        as :class:`pySecDec.algebra.Pow` with ``exponent = -1``
        and the ``base`` of type :class:`pySecDec.algebra.Polynomial`.

    :param indices:
        arbitrarily many integers;
        The index/indices of the parameter(s) to partially integrate.
        :math:`j` in the formulae above.

    Return the pole part and the numerically integrable remainder
    as a list. That is the sum and the integrand of equation (21)
    in arXiv:0803.4177v2 [Hei08]_.
    Each returned list element has the same structure as the input
    `polyprod`.

    '''
    new_products = [polyprod]
    for index in indices:
        old_products = new_products
        new_products = []
        for polyprod in old_products:
            new_products.extend( _integrate_pole_part_single_index(polyprod, index) )
    return new_products

# -------------------------------------------- integration by parts --------------------------------------------

def _integrate_by_parts_single_index(polyprod, power_goal, index):
    monomial_product = polyprod.factors[0] # e.g. (z1**2)**1-eps0-eps1  *  (z1*z2)**3-4*eps1
    regulator_poles = polyprod.factors[1]
    cal_I = Product(*polyprod.factors[2:], copy=False)
    polysymbols = monomial_product.symbols

    # extract the overall power from the `monomial_product`
    exponent_constant_term = 0
    full_exponent = 0
    for monomial_factor in monomial_product.factors:
        if monomial_factor.exponent.has_constant_term():
            this_factor_exponent_constant_term = np.sum(monomial_factor.exponent.coeffs[np.where((monomial_factor.exponent.expolist == 0).all(axis=1))])
            exponent_constant_term += this_factor_exponent_constant_term * monomial_factor.expolist[0,index]
        # else: e.g. (z1**2)**-eps0-eps1 --> no contribution to the regulator-free term
            # exponent_constant_term += 0
        full_exponent += monomial_factor.exponent * monomial_factor.expolist[0,index]

    if exponent_constant_term < 0:
        regulator_part = full_exponent - exponent_constant_term
        # the next line conceptually implements "assert regulator_part != 0"
        if (regulator_part.coeffs == 0).all():
            raise ValueError('"1/0" detected. An additional regulator is needed to calculate this integral.')

    def increase_monomial_power_by_one(monomial_product):
        # modify one of the factors if there is one raised to power ``1``
        have_suitable_factor = False
        for factor_index,factor in enumerate(monomial_product.factors):
            if (factor.exponent.expolist == 0).all() and (factor.exponent.coeffs == 1).all():
                have_suitable_factor = True
                break

        if have_suitable_factor:
            monomial_product = monomial_product.copy()
            factor = monomial_product.factors[factor_index]
            factor.expolist[0,index] += 1
            return monomial_product

        # fallback solution: add a new factor
        expolist = np.zeros((1, monomial_product.number_of_variables), dtype=int)
        expolist[:,index] = 1
        exponent = Polynomial(np.zeros_like(expolist,dtype=int), np.array([1]), polysymbols, copy=False)
        new_factor = ExponentiatedPolynomial(expolist, np.array([1]), exponent, polysymbols, copy=False)
        return Product(new_factor, *monomial_product.factors)

    # stop if `power_goal` is reached
    if exponent_constant_term >= power_goal:
        return [polyprod]

    # the factor ``1/(a+1-b*eps1-c*eps2-...)``
    regulator_denominator = full_exponent + 1
    new_pole_part = Pow(regulator_poles.base * regulator_denominator, regulator_poles.exponent, copy=False)

    # construct the term without integral
    term_without_integral_monomial_part = monomial_product.replace(index,1)
    term_without_integral_cal_I = cal_I.replace(index,1)
    term_without_integral = Product(term_without_integral_monomial_part, new_pole_part.copy(), term_without_integral_cal_I.simplify(), copy=False)

    # construct the term to be integrated
    term_with_integral_monomial_part = increase_monomial_power_by_one(monomial_product)
    term_with_integral_cal_I = -cal_I.derive(index)
    term_with_integral = Product(term_with_integral_monomial_part, new_pole_part.copy(), term_with_integral_cal_I.simplify(), copy=False)

    return [term_without_integral] + _integrate_by_parts_single_index(term_with_integral, power_goal, index)

def integrate_by_parts(polyprod, power_goals, indices):
    r'''
    Repeatedly apply integration by parts,

    .. math::
        \int_0^1
        {
            dt_j t_j^{(a_j - b_j \epsilon_1 - c \epsilon_2 + ...)}
            \mathcal{I} (t_j,\{t_{i \neq j}\}, \epsilon_1, \epsilon_2, ...)
        }
        =
        \frac{1}{a_j + 1 - b_j \epsilon_1 - c \epsilon_2 - ...}
        \left(
        \mathcal{I} (1,\{t_{i \neq j}\}, \epsilon_1, \epsilon_2, ...)
        -
        \int_0^1
        {
            dt_j t_j^{(a_j + 1 - b_j \epsilon_1 - c \epsilon_2 + ...)}
            \mathcal{I}^{\prime} (t_j,\{t_{i \neq j}\}, \epsilon_1, \epsilon_2, ...)
        }
        \right)

    , where :math:`\mathcal{I}^{\prime}` denotes the derivative
    of :math:`\mathcal{I}` with respect to :math:`t_j`. The iteration
    stops, when :math:`a_j>=` `power_goal_j`.

    .. seealso::
        This function provides an alternative to :func:`.integrate_pole_part`.

    :param polyprod:
        :class:`.algebra.Product` of the
        form ``<product of <monomial>**(a_j + ...)> *
        <regulator poles of cal_I> * <cal_I>``;
        The input product as decribed above.
        The <product of <monomial>**(a_j + ...)> should be
        a :class:`pySecDec.algebra.Product` of
        <monomial>**(a_j + ...).
        as described below.
        The <monomial>**(a_j + ...) should be an
        :class:`pySecDec.algebra.ExponentiatedPolynomial`
        with ``exponent`` being a :class:`.Polynomial` of the
        regulators :math:`\epsilon_1, \epsilon_2, ...`. Although
        no dependence on the Feynman parameters is expected
        in the ``exponent``, the polynomial variables should
        be the Feynman parameters and the regulators.
        The constant term of the exponent should be numerical.
        The polynomial variables of ``monomial`` and the other
        factors (interpreted as :math:`\mathcal{I}`) are interpreted
        as the Feynman parameters and the epsilon regulators.
        Make sure that the last factor (``<cal_I>``) is defined
        and finite for :math:`\epsilon = 0`. All poles for
        :math:`\epsilon \rightarrow 0` should be made explicit
        by putting them into ``<regulator poles of cal_I>``
        as :class:`pySecDec.algebra.Pow` with ``exponent = -1``
        and the ``base`` of type :class:`pySecDec.algebra.Polynomial`.

    :param power_goals:
        number or iterable of numbers, e.g. float, integer, ...;
        The stopping criterion for the iteration.

    :param indices:
        iterable of integers;
        The index/indices of the parameter(s) to partially integrate.
        :math:`j` in the formulae above.

    Return the pole part and the numerically integrable remainder
    as a list.
    Each returned list element has the same structure as the input
    `polyprod`.

    '''
    if not isinstance(indices,list):
        indices = list(indices)
    if np.iterable(power_goals):
        if not isinstance(power_goals,list):
            power_goals = list(power_goals)
        assert len(power_goals) == len(indices), 'The number of `power_goals` (%i) must equal the number of indices (%i).' % (len(power_goals), len(indices))
    else:
        power_goals = repeat(power_goals)
    new_products = [polyprod]
    for power_goal,index in zip(power_goals,indices):
        old_products = new_products
        new_products = []
        for polyprod in old_products:
            new_products.extend( _integrate_by_parts_single_index(polyprod, power_goal, index) )
    return new_products
