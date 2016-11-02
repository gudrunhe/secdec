"""Unit tests for the subtraction routines"""

from .subtraction import *
from .algebra import *
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
                    expected_pole_part = sp.sympify('''
                                                         (  1/(-2 + 0 + 1 - eps0 - 3*eps1) * (A + B*x0 + C*x0**2*x1 + D*x1) +
                                                            1/(-2 + 1 + 1 - eps0 - 3*eps1) * (B + 2*C*x0*x1)  ) * (x1**2)**(-2 - eps0 - 3*eps1)
                                                    ''').subs('x0',0)

                    expected_numerical_integrand = (sp.sympify( str(self.cal_I) ) - sp.sympify('A + D*x1 + x0*B') ) * sp.sympify( str(self.exponentiated_monomial1) )

                    if i == 0:
                        for summand in I_j_pole_part.summands:
                            np.testing.assert_array_equal(summand.factors[0].factors[0].expolist, [(0,2,0,0)])
                            np.testing.assert_array_equal(summand.factors[0].factors[0].coeffs, [1])
                            self.assertEqual(  sp.sympify( str(summand.factors[0].factors[0].exponent) )  -  sp.sympify('-2 - eps0 - 3*eps1'),  0  )
                    elif i == 1:
                        for summand in I_j_pole_part.summands:
                            np.testing.assert_array_equal(summand.factors[0].factors[0].expolist, [(0,2,0,0)])
                            np.testing.assert_array_equal(summand.factors[0].factors[1].expolist, [(0,4,0,0)])
                            np.testing.assert_array_equal(summand.factors[0].factors[0].coeffs, [1])
                            np.testing.assert_array_equal(summand.factors[0].factors[1].coeffs, [1])
                            self.assertEqual(  sp.sympify( str(summand.factors[0].factors[0].exponent) )  -  sp.sympify(-2),  0  )
                            self.assertEqual(  sp.sympify( str(summand.factors[0].factors[1].exponent) )  -  sp.sympify('- eps0 - 3*eps1')/2,  0  )
                    else:
                        raise IndexError('`i` should only run over `(0,1)`!')

                elif j == 1:
                    expected_pole_part = sp.sympify('''
                                                         (  1/( 2*(-2 - eps0 - 3*eps1) + 0 + 1 ) * (A + B*x0 + C*x0**2*x1 + D*x1) +
                                                            1/( 2*(-2 - eps0 - 3*eps1) + 1 + 1 ) * (C*x0**2 + D) +
                                                            0  ) * x0**(-2 - eps0 - 3*eps1)
                                                    ''').subs('x1',0)

                    expected_numerical_integrand = (sp.sympify( str(self.cal_I) ) - sp.sympify('A + B*x0 + C*x0**2*x1 + D*x1')) * sp.sympify( str(self.exponentiated_monomial1) )

                else:
                    raise IndexError('`j` should only run over `(0,1)`!')

                self.assertEqual( sp.powdenest((sp.sympify(str(I_j_pole_part)) - expected_pole_part), force=True).simplify() , 0)
                for summand in I_j_pole_part.summands:
                    self.assertEqual(type(summand), Product)
                    self.assertEqual(type(summand.factors[0]), Product)
                    for factor in summand.factors[0].factors:
                        self.assertEqual(type(factor), ExponentiatedPolynomial)

                self.assertEqual(type(I_j_numerically_integrable_part), Product)
                self.assertEqual(type(I_j_numerically_integrable_part.factors[0]), Product)
                for factor in I_j_numerically_integrable_part.factors[0].factors:
                    self.assertEqual(type(factor), ExponentiatedPolynomial)
                should_be_zero = (sp.sympify(I_j_numerically_integrable_part) - expected_numerical_integrand).simplify()
                if i == 1:
                    # need some simplifications that are not generally true for all complex numbers
                    # see http://docs.sympy.org/dev/tutorial/simplification.html#powers for a discussion
                    should_be_zero = sp.expand_power_base(should_be_zero, force=True)
                    should_be_zero = sp.powdenest(should_be_zero, force=True)
                self.assertEqual( should_be_zero , 0)

    def test_integrate_multiple_pole_parts(self):
        I_j_before = Product(self.monomial_product1,self.regulator_poles,self.cal_I)
        I_j_after = Sum(*integrate_pole_part(I_j_before,0,1))
        # expected_after_0 = sp.sympify('''
        #                                    1/(-2 + 0 + 1 - eps0 - 3*eps1) * (A + D*x1) * (x1)**(-4 - 2*eps0 - 6*eps1) +
        #                                    1/(-2 + 1 + 1 - eps0 - 3*eps1) * (B) * (x1)**(-4 - 2*eps0 - 6*eps1) +
        #                                    C*x0**2*x0**(-eps0 - 3*eps1 - 2) * (x1)**(-2*eps0 - 6*eps1 - 3)
        #                               ''')
        expected_after_0_1 = sp.sympify('''
                                             1/(-2 + 0 + 1 - eps0 - 3*eps1) / (-4+1-2*eps0-6*eps1) * (A) +
                                             1/(-2 + 0 + 1 - eps0 - 3*eps1) / (-4+2-2*eps0-6*eps1)* (D) +
                                             (  1/(-2 + 0 + 1 - eps0 - 3*eps1) * (A + D*x1)
                                                - 1/(-2 + 0 + 1 - eps0 - 3*eps1) * (A)
                                                - 1/(-2 + 0 + 1 - eps0 - 3*eps1) * (D) * x1
                                             ) * (x1)**(-4 - 2*eps0 - 6*eps1) +

                                             1/(-2 + 1 + 1 - eps0 - 3*eps1) / (-4+1-2*eps0-6*eps1) * (B) +

                                             C*x0**2*x0**(-eps0 - 3*eps1 - 2) / (-2*eps0 - 6*eps1 - 3+1)
                                        ''')
        self.assertEqual( (sp.sympify(str(I_j_after)) - expected_after_0_1).simplify() , 0)

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
