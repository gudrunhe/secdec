from .subtraction import *
from .algebra import *
from .misc import sympify_expression
import sympy as sp
import unittest
from nose.plugins.attrib import attr

#@attr('active')
class TestIntegratePolePart(unittest.TestCase):
    def setUp(self):
        self.Feynman_parameter_symbols = ['x0','x1']
        self.regulator_symbols = ['eps0','eps1']
        all_symbols = self.Feynman_parameter_symbols+self.regulator_symbols

        self.cal_I = Polynomial([(0,0,0,0),(1,0,0,0),(2,1,0,0),(0,1,0,0)],['A','B','C','D'], polysymbols=all_symbols)
        self.regulator_poles = Pow(Polynomial.from_expression('1', all_symbols), Expression(-1, all_symbols)) # initializer, no poles yet

        self.monomial_exponent1 = Polynomial.from_expression('-2 - eps0 - 3*eps1',self.Feynman_parameter_symbols+self.regulator_symbols)
        self.monomial_exponent2 = Polynomial.from_expression(-2,self.Feynman_parameter_symbols+self.regulator_symbols)
        self.monomial_exponent3 = Polynomial.from_expression('- 1/2*eps0 - 3/2*eps1',self.Feynman_parameter_symbols+self.regulator_symbols)
        self.exponentiated_monomial1 = ExponentiatedPolynomial([(1,2,0,0)],[1], exponent=self.monomial_exponent1, polysymbols=self.Feynman_parameter_symbols+self.regulator_symbols)
        self.exponentiated_monomial2 = ExponentiatedPolynomial([(1,2,0,0)],[1], exponent=self.monomial_exponent2, polysymbols=self.Feynman_parameter_symbols+self.regulator_symbols)
        self.exponentiated_monomial3 = ExponentiatedPolynomial([(2,4,0,0)],[1], exponent=self.monomial_exponent3, polysymbols=self.Feynman_parameter_symbols+self.regulator_symbols)

        self.monomial_product1 = Product(self.exponentiated_monomial1)
        self.monomial_product2 = Product(self.exponentiated_monomial2, self.exponentiated_monomial3)

    def test_integrate_pole_part(self):
        for j in (0,1):
            for i, I_j_before in enumerate([Product(self.monomial_product1,self.regulator_poles,self.cal_I),
                                            Product(self.monomial_product2,self.regulator_poles,self.cal_I)]):
                I_j_after = integrate_pole_part(I_j_before,j)
                I_j_pole_part = Sum(*I_j_after[:-1])
                I_j_numerically_integrable_part = I_j_after[-1]

                if j == 0:
                    expected_pole_part = sympify_expression('''
                                                         (  1/(-2 + 0 + 1 - eps0 - 3*eps1) * (A + B*x0 + C*x0**2*x1 + D*x1) +
                                                            1/(-2 + 1 + 1 - eps0 - 3*eps1) * (B + 2*C*x0*x1)  ) * (x1**2)**(-2 - eps0 - 3*eps1)
                                                    ''').subs('x0',0)

                    expected_numerical_integrand = (sympify_expression( str(self.cal_I) ) - sympify_expression('A + D*x1 + x0*B') ) * sympify_expression( str(self.exponentiated_monomial1) )

                    if i == 0:
                        for summand in I_j_pole_part.summands:
                            np.testing.assert_array_equal(summand.factors[0].factors[0].expolist, [(0,2,0,0)])
                            np.testing.assert_array_equal(summand.factors[0].factors[0].coeffs, [1])
                            self.assertEqual(  sympify_expression( str(summand.factors[0].factors[0].exponent) )  -  sympify_expression('-2 - eps0 - 3*eps1'),  0  )
                    elif i == 1:
                        for summand in I_j_pole_part.summands:
                            np.testing.assert_array_equal(summand.factors[0].factors[0].expolist, [(0,2,0,0)])
                            np.testing.assert_array_equal(summand.factors[0].factors[1].expolist, [(0,4,0,0)])
                            np.testing.assert_array_equal(summand.factors[0].factors[0].coeffs, [1])
                            np.testing.assert_array_equal(summand.factors[0].factors[1].coeffs, [1])
                            self.assertEqual(  sympify_expression( str(summand.factors[0].factors[0].exponent) )  -  sympify_expression(-2),  0  )
                            self.assertEqual(  sympify_expression( str(summand.factors[0].factors[1].exponent) )  -  sympify_expression('- eps0 - 3*eps1')/2,  0  )
                    else:
                        raise IndexError('`i` should only run over `(0,1)`!')

                elif j == 1:
                    expected_pole_part = sympify_expression('''
                                                         (  1/( 2*(-2 - eps0 - 3*eps1) + 0 + 1 ) * (A + B*x0 + C*x0**2*x1 + D*x1) +
                                                            1/( 2*(-2 - eps0 - 3*eps1) + 1 + 1 ) * (C*x0**2 + D) +
                                                            0  ) * x0**(-2 - eps0 - 3*eps1)
                                                    ''').subs('x1',0)

                    expected_numerical_integrand = (sympify_expression( str(self.cal_I) ) - sympify_expression('A + B*x0 + C*x0**2*x1 + D*x1')) * sympify_expression( str(self.exponentiated_monomial1) )

                else:
                    raise IndexError('`j` should only run over `(0,1)`!')

                self.assertEqual( sp.powdenest((sympify_expression(str(I_j_pole_part)) - expected_pole_part), force=True).simplify() , 0)
                for summand in I_j_pole_part.summands:
                    self.assertEqual(type(summand), Product)
                    self.assertEqual(type(summand.factors[0]), Product)
                    for factor in summand.factors[0].factors:
                        self.assertEqual(type(factor), ExponentiatedPolynomial)

                self.assertEqual(type(I_j_numerically_integrable_part), Product)
                self.assertEqual(type(I_j_numerically_integrable_part.factors[0]), Product)
                for factor in I_j_numerically_integrable_part.factors[0].factors:
                    self.assertEqual(type(factor), ExponentiatedPolynomial)
                should_be_zero = (sympify_expression(I_j_numerically_integrable_part) - expected_numerical_integrand).simplify()
                if i == 1:
                    # need some simplifications that are not generally true for all complex numbers
                    # see http://docs.sympy.org/dev/tutorial/simplification.html#powers for a discussion
                    should_be_zero = sp.expand_power_base(should_be_zero, force=True)
                    should_be_zero = sp.powdenest(should_be_zero, force=True)
                self.assertEqual( should_be_zero , 0)

    def test_integrate_multiple_pole_parts(self):
        I_j_before = Product(self.monomial_product1,self.regulator_poles,self.cal_I)
        I_j_after = Sum(*integrate_pole_part(I_j_before,0,1))
        # expected_after_0 = sympify_expression('''
        #                                    1/(-2 + 0 + 1 - eps0 - 3*eps1) * (A + D*x1) * (x1)**(-4 - 2*eps0 - 6*eps1) +
        #                                    1/(-2 + 1 + 1 - eps0 - 3*eps1) * (B) * (x1)**(-4 - 2*eps0 - 6*eps1) +
        #                                    C*x0**2*x0**(-eps0 - 3*eps1 - 2) * (x1)**(-2*eps0 - 6*eps1 - 3)
        #                               ''')
        expected_after_0_1 = sympify_expression('''
                                             1/(-2 + 0 + 1 - eps0 - 3*eps1) / (-4+1-2*eps0-6*eps1) * (A) +
                                             1/(-2 + 0 + 1 - eps0 - 3*eps1) / (-4+2-2*eps0-6*eps1)* (D) +
                                             (  1/(-2 + 0 + 1 - eps0 - 3*eps1) * (A + D*x1)
                                                - 1/(-2 + 0 + 1 - eps0 - 3*eps1) * (A)
                                                - 1/(-2 + 0 + 1 - eps0 - 3*eps1) * (D) * x1
                                             ) * (x1)**(-4 - 2*eps0 - 6*eps1) +

                                             1/(-2 + 1 + 1 - eps0 - 3*eps1) / (-4+1-2*eps0-6*eps1) * (B) +

                                             C*x0**2*x0**(-eps0 - 3*eps1 - 2) / (-2*eps0 - 6*eps1 - 3+1)
                                        ''')
        self.assertEqual( (sympify_expression(str(I_j_after)) - expected_after_0_1).simplify() , 0)

    def test_subsequent_call(self):
        # should be able to walk through all variables
        I = Product(self.monomial_product1,self.regulator_poles,self.cal_I)
        I_0 = integrate_pole_part(I,0)
        for i,term in enumerate(I_0):
            integrate_pole_part(term,1)

    #@attr('active')
    def test_catch_one_over_zero(self):
        minus_one = Polynomial([[0]], [-1])
        monomial = ExponentiatedPolynomial([[1]], [1], exponent=minus_one) # "x0**(-1)" without regulator --> leads to "1/0" in subtraction
        pole_part_initializer = Pow(Polynomial([[0]],[1]), exponent=minus_one)
        cal_I = Function('cal_I', Polynomial([[1]],[1])) # "cal_I(x0)"
        prod = Product(Product(monomial), pole_part_initializer, cal_I)
        self.assertRaisesRegexp(ValueError, '1/0', integrate_pole_part, prod, 0)

    #@attr('active')
    def test_pole_structure(self):
        self.assertEqual( pole_structure(self.monomial_product1,0,1) , [-2,-4] )
        self.assertEqual( pole_structure(self.monomial_product2,0,1) , [-2,-4] )

#@attr('active')
class TestIntegrateByParts(unittest.TestCase):
    def setUp(self):
        self.symbols = ['x0', 'x1', 'x2', 'eps1', 'eps2']
        a, b = sp.symbols('a b')
        eps1 = Polynomial.from_expression('eps1', self.symbols)
        eps2 = Polynomial.from_expression('eps2', self.symbols)
        polynomial_one = Polynomial(np.zeros([1,len(self.symbols)], dtype=int), np.array([1]), self.symbols, copy=False)

        def xi_to_the(index, power):
            expovec = [0] * len(self.symbols)
            expovec[index] = 1
            return ExponentiatedPolynomial([expovec], [1], polysymbols=self.symbols, exponent=power + a*eps1 + b*eps2)

        # construct `ibp_input` as:
        # ``x0**(-1 + a*eps1 + b*eps2) * x1**(-2 + a*eps1 + b*eps2) * x2**(2 + a*eps1 + b*eps2) * cal_I(x0,x1,x2,eps1,eps2)``
        monomials = Product( xi_to_the(0,-1), xi_to_the(1,-2), xi_to_the(2,2) )
        pole_part_initializer = Pow(polynomial_one, -polynomial_one)
        cal_I = Function(   'cal_I' , *(  Polynomial.from_expression(symbol, self.symbols) for symbol in self.symbols  )   )

        self.ibp_input = Product(monomials, pole_part_initializer, cal_I)

    def check_terms(self, computed, expected):
        self.assertEqual( len(computed) , len(expected) )
        for index,(term,target_term) in enumerate(zip(computed, expected)):
            print('index: %i' % index)
            sympified_term = sympify_expression(term)
            sympified_target_term = sympify_expression(target_term)
            difference = (sympified_term - sympified_target_term).simplify()
            self.assertEqual(difference,0)

    #@attr('active')
    def test_error_messages(self):
        self.assertRaisesRegexp(AssertionError, 'number of .*power_goals.* \(1\).* (match|equal) .*number of .*indices.* \(2\)', integrate_by_parts, self.ibp_input, [-1], [0,2])

    #@attr('active')
    def test_power_goal_minus1(self):
        terms_after_ibp = integrate_by_parts(self.ibp_input, -1, [0,1,2])
        target_terms_after_ibp = \
        [
            'x0**(-1 + a*eps1 + b*eps2) * x2**(2 + a*eps1 + b*eps2) * 1/(-1 + a*eps1 + b*eps2) * cal_I(x0,1,x2,eps1,eps2)',
            '- x0**(-1 + a*eps1 + b*eps2) * x1**(-1 + a*eps1 + b*eps2) * x2**(2 + a*eps1 + b*eps2) * 1/(-1 + a*eps1 + b*eps2) * dcal_Id1(x0,x1,x2,eps1,eps2)'
        ]
        self.check_terms(terms_after_ibp, target_terms_after_ibp)

    #@attr('active')
    def test_power_goal_0(self):
        terms_after_ibp_input_together = integrate_by_parts(self.ibp_input, 0, (0,1,2))
        terms_after_ibp_input_separate = integrate_by_parts(self.ibp_input, (0,0,0), (0,1,2))
        target_terms_after_ibp = \
        [
            'x2**(2 + a*eps1 + b*eps2) * 1/(a*eps1 + b*eps2) * 1/(-1 + a*eps1 + b*eps2) * cal_I(1,1,x2,eps1,eps2)',
            '- x2**(2 + a*eps1 + b*eps2) * 1/(a*eps1 + b*eps2)**2 * 1/(-1 + a*eps1 + b*eps2) * dcal_Id1(1,1,x2,eps1,eps2)',
            'x1**(a*eps1 + b*eps2) * x2**(2 + a*eps1 + b*eps2) * 1/(a*eps1 + b*eps2)**2 * 1/(-1 + a*eps1 + b*eps2) * ddcal_Id1d1(1,x1,x2,eps1,eps2)',
            'x2**(2 + a*eps1 + b*eps2) * 1/(a*eps1 + b*eps2) * 1/(-1 + a*eps1 + b*eps2) * (- x0**(a*eps1 + b*eps2) * dcal_Id0(x0,1,x2,eps1,eps2))',
            '- x2**(2 + a*eps1 + b*eps2) * 1/(a*eps1 + b*eps2)**2 * 1/(-1 + a*eps1 + b*eps2) * (- x0**(a*eps1 + b*eps2) * ddcal_Id0d1(x0,1,x2,eps1,eps2))',
            'x1**(a*eps1 + b*eps2) * x2**(2 + a*eps1 + b*eps2) * 1/(a*eps1 + b*eps2)**2 * 1/(-1 + a*eps1 + b*eps2) * (- x0**(a*eps1 + b*eps2) * dddcal_Id0d1d1(x0,x1,x2,eps1,eps2))'
        ]
        self.check_terms(terms_after_ibp_input_together, target_terms_after_ibp)
        self.check_terms(terms_after_ibp_input_separate, target_terms_after_ibp)

    #@attr('active')
    def test_power_goals_0_minus1_1(self):
        terms_after_ibp = integrate_by_parts(self.ibp_input, [0,-1,1], [0,1,2])
        target_terms_after_ibp = \
        [
            'x2**(a*eps1 + b*eps2 + 2)*cal_I(1, 1, x2, eps1, eps2)/(a**2*eps1**2 + 2*a*b*eps1*eps2 - a*eps1 + b**2*eps2**2 - b*eps2)',
            '-x1*x1**(a*eps1 + b*eps2 - 2)*x2**(a*eps1 + b*eps2 + 2)*dcal_Id1(1, x1, x2, eps1, eps2)/(a**2*eps1**2 + 2*a*b*eps1*eps2 - a*eps1 + b**2*eps2**2 - b*eps2)',
            '-x0*x0**(a*eps1 + b*eps2 - 1)*x2**(a*eps1 + b*eps2 + 2)*dcal_Id0(x0, 1, x2, eps1, eps2)/(a**2*eps1**2 + 2*a*b*eps1*eps2 - a*eps1 + b**2*eps2**2 - b*eps2)',
            'x0*x0**(a*eps1 + b*eps2 - 1)*x1*x1**(a*eps1 + b*eps2 - 2)*x2**(a*eps1 + b*eps2 + 2)*ddcal_Id0d1(x0, x1, x2, eps1, eps2)/(a**2*eps1**2 + 2*a*b*eps1*eps2 - a*eps1 + b**2*eps2**2 - b*eps2)'
        ]
        self.check_terms(terms_after_ibp, target_terms_after_ibp)

    #@attr('active')
    def test_select_index(self):
        terms_after_ibp = integrate_by_parts(self.ibp_input, 0, [0])
        target_terms_after_ibp = \
        [
            'x1**(-2 + a*eps1 + b*eps2) * x2**(2 + a*eps1 + b*eps2) * 1/(a*eps1 + b*eps2) * cal_I(1,x1,x2,eps1,eps2)',
            'x1**(-2 + a*eps1 + b*eps2) * x2**(2 + a*eps1 + b*eps2) * 1/(a*eps1 + b*eps2) * (- x0**(a*eps1 + b*eps2) * dcal_Id0(x0,x1,x2,eps1,eps2))'
        ]
        self.check_terms(terms_after_ibp, target_terms_after_ibp)
