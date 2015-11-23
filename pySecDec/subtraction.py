"""
Routines to isolate the divergencies in an :math:`\epsilon`
expansion

"""

from .polynomial import Polynomial, ExponentiatedPolynomial, PolynomialSum, PolynomialProduct, replace
import numpy as np
import sympy as sp

def integrate_pole_part(polyprod, index):
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
        :class:`.polynomial.PolynomialProduct` of the
        form ``<monomial>**(a_j + ...) * <regulator poles of cal_I>
        * <cal_I>``;
        The input product as decribed above.
        The monomial should be an
        :class:`pySecDec.polynomial.ExponentiatedPolynomial`
        with ``exponent`` being a :class:`Polynomial` of the
        regulators :math:`\epsilon_1, \epsilon_2, ...`.
        The constant term of the exponent should be numerical.
        The polynomial variables of ``monomial`` and the other
        factors (interpreted as :math:`\cal{I}`) are interpreted
        as the Feynman parameters and the epsilon regulators.

    :param index:
        integer;
        The index of the parameter to partially integrate.
        :math:`j` in the formulae above.

    Return the pole part and the numerically integrable remainder
    as a list. That is the sum and the integrand of equation (21)
    in arXiv:0803.4177v2.
    Each returned list element has the same structure as the input
    `polyprod`.

    '''
    monomial = polyprod.factors[0]
    monomial_FeynmanJ_set_to_one = replace(monomial,index,1)
    regulator_poles = polyprod.factors[1]
    cal_I = PolynomialProduct(*polyprod.factors[2:])
    polysymbols = monomial.polysymbols
    # arXiv:0803.4177v2: exponent_constant_term = a_j
    if monomial.exponent.has_constant_term():
        exponent_constant_term  = np.sum(monomial.exponent.coeffs[np.where((monomial.exponent.expolist == 0).all(axis=1))])
        exponent_constant_term *= monomial.expolist[0,index]
    else:
        exponent_constant_term = 0
    number_of_regulators = monomial.exponent.number_of_variables
    number_of_Feynman_parameters = monomial.number_of_variables - number_of_regulators
    assert number_of_Feynman_parameters > 0, 'unexpected input'
    assert type(monomial.exponent) is Polynomial , 'unexpected input'

    # convert the exponent to a polynomial that also has slots for the Feynman parameters
    monomial_exponent_expolist_with_Feynman_parameters = \
    np.hstack([ np.zeros([len(monomial.exponent.expolist),number_of_Feynman_parameters], dtype=int) , monomial.exponent.expolist ])
    monomial_exponent_with_Feynman_parameters = Polynomial(monomial_exponent_expolist_with_Feynman_parameters, monomial.exponent.coeffs, polysymbols=polysymbols)

    def make_FeynmanIndex_to_power(power):
        '''
        Return a :class:`Polynomial` representing
        "Feynman_index**power/power!"

        '''
        expolist = [0]*monomial.number_of_variables
        expolist[index] = power
        return Polynomial([expolist], ['1/%s!'%power])


    if exponent_constant_term > -1:
        # no subtraction needed, the input `polyprod` is numerically integrable
        return [polyprod]


    output_summands = []
    minus_cal_I_expansion_summands = []
    # exponent_constant_term <= -1 < 0      =>     abs(exponent_constant_term) = -exponent_constant_term
    # use symbol names as in arXiv:0803.4177v2

    # construction of the pole part
    for p in range(int(-exponent_constant_term)):
        # No dependency on Feynman parameter with index `index` --> has been integrated out analytically
        #   --> Replace that Feynman parameter by `1` in the monomial factor.
        current_factors = [monomial_FeynmanJ_set_to_one]

        # renew `p_factorial` on the fly
        p_factorial = 1 if p == 0 else p_factorial * p

        # renew derivative of cal_I on the fly
        derivative_cal_I = cal_I if p == 0 else derivative_cal_I.derive(index)
        derivative_cal_I_Feynmanj_set_to_zero = replace(derivative_cal_I, index, 0)

        # arXiv0803.4177v2: 1/( (a_j + p + 1 - b_j * eps) * factorial(p) )
        new_potential_pole_denominator = (monomial.expolist[0,index] * monomial_exponent_with_Feynman_parameters + p + 1) * (p_factorial)
        new_potential_pole = ExponentiatedPolynomial(new_potential_pole_denominator.expolist, new_potential_pole_denominator.coeffs, exponent=-1, polysymbols=monomial.polysymbols)
        # put this factor into the pole part only if a_j + p is zero
        if exponent_constant_term + p + 1 == 0:
            current_regulator_poles = PolynomialProduct(new_potential_pole, regulator_poles).simplify()
            current_factors.append(current_regulator_poles)
            current_factors.append(derivative_cal_I_Feynmanj_set_to_zero)
        # otherwise it does not lead to additional regulator poles and can become part of <derivative_cal_I>
        else:
            # poles in current term: none --> poles are just the old ones
            current_factors.append(regulator_poles)

            # put `new_potential_pole` (which is not a pole in this case) in the last factor
            last_factor = PolynomialProduct(derivative_cal_I_Feynmanj_set_to_zero, new_potential_pole)
            current_factors.append(last_factor.simplify())

        output_summands.append(PolynomialProduct(*current_factors))
        minus_cal_I_expansion_summands.append(PolynomialProduct(-make_FeynmanIndex_to_power(p), derivative_cal_I_Feynmanj_set_to_zero ))

    integrable_part = PolynomialProduct(   monomial, regulator_poles, PolynomialSum(cal_I, *minus_cal_I_expansion_summands).simplify()   )

    output_summands.append(integrable_part)

    return output_summands
