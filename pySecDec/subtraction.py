"""
Subtraction
-----------

Routines to isolate the divergencies in an :math:`\epsilon`
expansion.

"""

from .algebra import Polynomial, Sum, Product, Pow
import numpy as np
import sympy as sp

_sympy_one = sp.sympify(1)

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
            raise ValueError('"1/0" detected.')

    # finite variance of Monte Carlo integral estimator only if ``exponent_constant_term > -0.5``
    # to be absolutely safe here, we eliminate integrable singularities; i.e. ``exponent_constant_term >= 0``?
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
            \cal{I} (t_j,\{t_{i \neq j}\}, \epsilon_1, \epsilon_2, ...)
        }

    into the form

    .. math::
        \sum_{p=0}^{|a|-1}
        {
            \frac{1}{a + p + 1 - b \epsilon_1 - c \epsilon_2 - ...}
            \frac{\cal{I}^{(p)} (0,\{t_{i \neq j}\}, \epsilon_1, \epsilon_2, ...)}{p!} +
            \int_0^1
            {
                dt_j t_j^{(a - b \epsilon_1 - c \epsilon_2 + ...)}
                R(t_j,\{t_{i \neq j}\}, \epsilon_1, \epsilon_2, ...)
            }
        }

    , where :math:`\cal{I}^{(p)}` is denotes the p-th derivative
    of :math:`\cal{I}` with respect to :math:`t_j`. The equations
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
        factors (interpreted as :math:`\cal{I}`) are interpreted
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
