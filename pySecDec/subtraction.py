"""
Routines to isolate the divergencies in an :math:`\epsilon`
expansion

"""

from .algebra import Polynomial, ExponentiatedPolynomial, Sum, Product
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
        assert type(monomial_factor.exponent) is Polynomial , 'unexpected input'
        if monomial_factor.exponent.has_constant_term():
            this_factor_exponent_constant_term = np.sum(monomial_factor.exponent.coeffs[np.where((monomial_factor.exponent.expolist == 0).all(axis=1))])
            exponent_constant_term += this_factor_exponent_constant_term * monomial_factor.expolist[0,index]
        # else: e.g. (z1**2)**-eps0-eps1 --> no contribution to ``a_j``
            # exponent_constant_term += 0
        full_exponent += monomial_factor.exponent * monomial_factor.expolist[0,index]

    # TODO: finite variance of Monte Carlo integral estimator only if exponent_constant_term >-0.5   -->   should we change this to >-0.5 or even >=0?
    if exponent_constant_term > -1:
        # no subtraction needed, the input `polyprod` is numerically integrable
        return [polyprod]

    def make_FeynmanIndex_to_power(power, factorial_of_power):
        '''
        Return a :class:`Polynomial` representing
        "Feynman_index**power/power!"

        '''
        expolist = [0]*monomial_product.number_of_variables
        expolist[index] = power
        return Polynomial([expolist], [_sympy_one/factorial_of_power], copy=False)

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
        new_potential_pole = ExponentiatedPolynomial(new_potential_pole_denominator.expolist, new_potential_pole_denominator.coeffs, exponent=-1, polysymbols=polysymbols, copy=False)
        # put this factor into the pole part only if a_j + p + 1 is zero
        if exponent_constant_term + p + 1 == 0:
            current_regulator_poles = Product(new_potential_pole, regulator_poles, copy=False)
            current_factors.append(current_regulator_poles)
            current_factors.append(derivative_cal_I_Feynmanj_set_to_zero)
        # otherwise it does not lead to additional regulator poles and can become part of <derivative_cal_I>
        else:
            # poles in current term: none --> poles are just the old ones
            current_factors.append(regulator_poles)

            # put `new_potential_pole` (which is not a pole in this case) in the last factor
            last_factor = Product(derivative_cal_I_Feynmanj_set_to_zero, new_potential_pole, copy=False)
            current_factors.append(last_factor.simplify())

        output_summands.append(Product(*current_factors, copy=False))
        minus_cal_I_expansion_summands.append(Product(-make_FeynmanIndex_to_power(p, p_factorial), derivative_cal_I_Feynmanj_set_to_zero, copy=False ))

    integrable_part = Product(   monomial_product, regulator_poles, Sum(cal_I, *minus_cal_I_expansion_summands, copy=False).simplify() , copy=False   )

    output_summands.append(integrable_part)

    return output_summands

# TODO: make this a private function --> should have a high-level function that takes a :class:`pySecDec.decomposition.Sector`.
# TODO: raise an error if the exponent depends on the Feynman parameters in that high-level function
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
        equation (19) to (21) as described in arXiv:0803.4177v2.


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
        with ``exponent`` being a :class:`Polynomial` of the
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
        as :class:`pySecDec.algebra.ExponentiatedPolynomial`
        with ``exponent = -1``.

    :param indices:
        arbitrarily many integers;
        The index/indices of the parameter(s) to partially integrate.
        :math:`j` in the formulae above.

    Return the pole part and the numerically integrable remainder
    as a list. That is the sum and the integrand of equation (21)
    in arXiv:0803.4177v2.
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


# TODO: implement something like this (including ``U`` as parameter):
#In [37]: def sum_gen(nunu):
#    for exponents, coeff in zip(nunu.expolist, nunu.coeffs):
#        exponents = list(exponents) + [0,0]
#        expr = Polynomial([exponents], [1], nunu.polysymbols + sp.sympify(['F','eps']))
#        yield Product( expr , make_expr(coeff, nunu.polysymbols + sp.sympify(['F','eps'])) )
#   ....:         
#
#In [38]: %time nun = Sum(*[item for item in sum_gen(li.numerator) ] )
#CPU times: user 366 ms, sys: 0 ns, total: 366 ms
#Wall time: 361 ms
#
