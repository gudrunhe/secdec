import pySecDec as psd
import sympy as sp
import numpy as np
from itertools import chain

# 2loop box with numerator
propagators = ['k1**2','(k1+p2)**2','(k1-p1)**2','(k1-k2)**2','(k2+p2)**2','(k2-p1)**2','(k2+p2+p3)**2']
numerator = 'k1(mu)*k1(mu) + 2*k1(mu)*p3(mu) + p3(mu)*p3(mu)' # (k1 + p3) ** 2
loop_momenta = ['k1','k2']
external_momenta = ['p1','p2','p3','p4']

replacement_rules = [
                        ('p1*p1', 0),
                        ('p2*p2', 0),
                        ('p3*p3', 0),
                        ('p4*p4', 0),
                        ('p1*p2', 's/2'),
                        ('p2*p3', 't/2'),
                        ('p1*p3', '-s/2-t/2')
                    ]

li = psd.loop_integral.LoopIntegral.from_propagators(propagators, loop_momenta, external_momenta,
                                                     numerator=numerator, replacement_rules=replacement_rules,
                                                     Feynman_parameters=['z%i'%i for i in range(1,7+1)],
                                                     Lorentz_indices=['mu'])

# the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
requested_order = 0

# highest pole from ``li.Gamma_factor`` and ``li.regulator ** li.regulator_power``
highest_prefactor_pole_order = - li.regulator_power
# ``li.Gamma_factor`` only leads to a pole if the argument for ``epsilon=0`` is less or equal to zero
if li.Gamma_factor.args[0].subs(li.regulator, 0) <= 0:
    highest_prefactor_pole_order += 1

# TODO: remove all sympy symbols from coeffs

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------

# hide ``U`` and ``F`` from the `polysymbols` of the `numerator`
# Having them as polysymbols can be useful if the part of numerator that does not depend on them shall be decomposed
# (that is sufficient for a complete decomposition of the numerator if ``U`` and ``F`` are decomposed as well)
numerator, numerator_F_U = psd.decomposition.hide(li.numerator, 2)

U = li.exponentiated_U
F = li.exponentiated_F

# we will not fully decompose the numerator, only ``F`` and ``U``
initial_sector = psd.decomposition.Sector([F, U], [numerator])

# decompose
# for the geometric/iterative method, change the following two calls accordingly
#primary_sectors = psd.decomposition.iterative.primary_decomposition(initial_sector)
#sector_generator = psd.decomposition.iterative.iterative_decomposition
#geomethod = False
primary_sectors = [psd.decomposition.geometric.Cheng_Wu(initial_sector)]
sector_generator = psd.decomposition.geometric.geometric_decomposition
geomethod = True


# subtraction and expansion
for primary_sector_index, primary_sector in enumerate(primary_sectors):
    if geomethod:
        sympy_params_without_regulator = li.Feynman_parameters[:-1]
    else:
        sympy_params_without_regulator = li.Feynman_parameters[:primary_sector_index] + li.Feynman_parameters[primary_sector_index+1:]
    sympy_params = sympy_params_without_regulator + [li.regulator]
    params_withoput_regulator = [psd.algebra.Polynomial.from_expression(p, sympy_params) for p in sympy_params_without_regulator]
    params = params_withoput_regulator + [psd.algebra.Polynomial.from_expression(li.regulator, sympy_params)]

    # define the symbolic functions ``F``, ``U``, and ``numerator``
    symbolic_F = psd.algebra.Function('F', *params_withoput_regulator)
    symbolic_U = psd.algebra.Function('U', *params_withoput_regulator)
    symbolic_numerator = psd.algebra.Function('numerator', *params)

    F_monomial, F_polynomial = primary_sector.cast[0].factors
    U_monomial, U_polynomial = primary_sector.cast[1].factors
    Jacobian = primary_sector.Jacobian

    # insert ``U`` and ``F`` as dummy functions into the numerator
    # note that ``U`` and ``F`` do not depend on the regulator
    numerator = psd.decomposition.unhide(primary_sector.other[0], numerator_F_U)
    numerator = numerator.replace(-1, sp.sympify('U')(*sympy_params_without_regulator), remove=True)
    numerator = numerator.replace(-1, sp.sympify('F')(*sympy_params_without_regulator), remove=True)

    # convert the coefficients and exponents to `_Expression`s
    hidden_variables = []
    for poly in [F_monomial, F_polynomial, U_monomial, U_polynomial, numerator, Jacobian]:
        new_coeffs = np.empty_like(poly.coeffs, dtype=object)
        for i, coeff in enumerate(poly.coeffs):
            new_coeffs[i] = psd.algebra.Expression(coeff, sympy_params)
        poly.coeffs = new_coeffs

        # append regulator to polysymbols
        poly.polysymbols.append(li.regulator)
        poly.number_of_variables += 1
        poly.expolist = np.hstack([poly.expolist, np.zeros([len(poly.expolist),1], dtype=int)])

        try:
            poly.exponent = psd.algebra.Expression(poly.exponent, sympy_params)
        except AttributeError:
            pass

        # hide regulator
        remainder, hidden = psd.decomposition.hide(poly, 1)
        hidden_variables.append(hidden)

    # apply changes to the primary_sector
    primary_sector = psd.decomposition.Sector([psd.algebra.Product(F_monomial, F_polynomial), psd.algebra.Product(U_monomial, U_polynomial)], [numerator], Jacobian)

    polynomial_one = psd.algebra.Polynomial(np.zeros([1,len(params)], dtype=int), np.array([1]), sympy_params, copy=False)
    pole_part_initializer = psd.algebra.Pow(polynomial_one, -polynomial_one) # this is just one packed into a suitable expression

    for sector_index, sector in enumerate(sector_generator(primary_sector)):
        F_monomial, F_polynomial = sector.cast[0].factors
        U_monomial, U_polynomial = sector.cast[1].factors
        unfactorized_numerator = sector.other[0]
        Jacobian = sector.Jacobian

        # unhide regulator
        for poly, hidden in zip([F_monomial, F_polynomial, U_monomial, U_polynomial, unfactorized_numerator, Jacobian], hidden_variables):
            psd.decomposition.unhide(poly, hidden)

        # insert the monomials of ``U`` and ``F`` into the numerator ``U/F --> xi**pow_i * U/F``
        # ``[:,:-1]`` is the part of the expolist that contains the Feynman parameter powers for `U_monomial`, `F_monomial`, and `unfactorized_numerator`
        # can do the replacement ``U/F --> xi**pow_i * U/F`` as follows:
        unfactorized_numerator.expolist[:,:-1] += np.einsum('i,k->ik', numerator_F_U.expolist[:,0], U_monomial.expolist[0,:-1]) # F
        unfactorized_numerator.expolist[:,:-1] += np.einsum('i,k->ik', numerator_F_U.expolist[:,1], F_monomial.expolist[0,:-1]) # U

        # factorize the numerator
        numerator = psd.algebra.Product(polynomial_one, unfactorized_numerator)
        psd.decomposition.refactorize(numerator)
        numerator_monomial, numerator_polynomial = numerator.factors

        # subtraction needs type `ExponentiatedPolynomial` in its monomial part
        Jacobian = psd.algebra.ExponentiatedPolynomial(sector.Jacobian.expolist, sector.Jacobian.coeffs, polysymbols=Jacobian.polysymbols, exponent=polynomial_one, copy=False)
        numerator_monomial = psd.algebra.ExponentiatedPolynomial(numerator_monomial.expolist, numerator_monomial.coeffs, polysymbols=numerator_monomial.polysymbols, exponent=polynomial_one, copy=False)

        # initialize the Product to be passed to `integrate_pole_part`
        monomials = psd.algebra.Product(Jacobian, F_monomial, U_monomial, numerator_monomial)
        #cal_I = psd.algebra.Product(F_polynomial, U_polynomial, numerator_polynomial) # explicitly inserted F and U # TODO: how to insert them in the numerator?
        #cal_I = psd.algebra.Product(psd.algebra.Pow(symbolic_F, F_polynomial.exponent), psd.algebra.Pow(symbolic_U, U_polynomial.exponent), numerator_polynomial) # symbolic F and U # TODO: how to find out which derivatives of them are needed
        cal_I = psd.algebra.Product(psd.algebra.Pow(symbolic_F, F_polynomial.exponent), psd.algebra.Pow(symbolic_U, U_polynomial.exponent), symbolic_numerator) # symbolic F, U and numerator # TODO: how to find out which derivatives of them are needed
        subtraction_initializer = psd.algebra.Product(monomials, pole_part_initializer, cal_I)

        subtracted = psd.subtraction.integrate_pole_part(subtraction_initializer, *range(subtraction_initializer.number_of_variables - 1))

        # expansion
        pole_parts = [s.factors[1].simplify() for s in subtracted]
        regular_parts = [psd.algebra.Product( *([s.factors[0]] + s.factors[2:]) ) for s in subtracted]

        for regular, singular in zip(regular_parts, pole_parts):
            # must expand every term to the requested order plus the highest pole it multiplies
            # We calculated the highest pole order of the prefactor (variable ``sector_epsilon_expansion_order``) above.
            # In addition, we have to take the poles of the current term into account.
            singular_expanded = psd.expansion.expand_singular(psd.algebra.Product(singular), -1, requested_order + highest_prefactor_pole_order)

            highest_pole_current_term = - singular_expanded.expolist[0,-1]
            sector_epsilon_expansion_order = requested_order + highest_prefactor_pole_order + highest_pole_current_term
            expansion_order = sector_epsilon_expansion_order
            regular_expanded = psd.expansion.expand_Taylor(regular, -1, sector_epsilon_expansion_order)

            # print singular_expanded * regular_expanded
            print 'expansion_order', expansion_order

            print 'sector_index', sector_index

        # TODO: which derivatives of U and F are needed? --> just calculate some for more realistic timings
        numerator_polynomial.derive(0).derive(1).derive(2)
        U_polynomial.derive(0).derive(1).derive(2)
        F_polynomial.derive(0).derive(1).derive(2)

        #if sector_index == 100:
        #    exit()

        # str_pole_p = str(pole_parts)

        #if '0))**' in str_pole_p:
        #    print
        #    print '******************'
        #    print '*****ALERT********'
        #    print '******************'
        #    print
        #    raise ValueError()
        # print str_pole_p

        #if sector_index == 10:
        #    exit()


