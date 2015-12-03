"""Unit tests for the subtraction routines"""

from .subtraction import *
from .algebra import *
import sympy as sp
import unittest

class TestIntegratePolePart(unittest.TestCase):
    def setUp(self):
        self.Feynman_parameter_symbols = ['x0','x1']
        self.regulator_symbols = ['eps0','eps1']

        self.cal_I = Polynomial([(0,0,0,0),(1,0,0,0),(2,1,0,0),(0,1,0,0)],['A','B','C','D'], polysymbols=self.Feynman_parameter_symbols+self.regulator_symbols)
        self.regulator_poles = Polynomial.from_expression('1', ['x0','x1','eps0','eps1']) # initializer, no poles yet

        self.monomial_exponent = Polynomial.from_expression('-2 - eps0 - 3*eps1',self.Feynman_parameter_symbols+self.regulator_symbols)
        self.exponentiated_monomial = ExponentiatedPolynomial([(1,2,0,0)],[1], exponent=self.monomial_exponent, polysymbols=self.Feynman_parameter_symbols+self.regulator_symbols)

    def test_integrate_pole_part(self):
        for j in (0,1):
            I_j_before = Product(self.exponentiated_monomial,self.regulator_poles,self.cal_I)
            I_j_after = integrate_pole_part(I_j_before,j)
            I_j_pole_part = PolynomialSum(*I_j_after[:-1])
            I_j_numerically_integrable_part = I_j_after[-1]

            if j == 0:
                expected_pole_part = sp.sympify('''
                                                     (  1/(-2 + 0 + 1 - eps0 - 3*eps1) * (A + B*x0 + C*x0**2*x1 + D*x1) +
                                                        1/(-2 + 1 + 1 - eps0 - 3*eps1) * (B + 2*C*x0*x1)  ) * (x1**2)**(-2 - eps0 - 3*eps1)
                                                ''').subs('x0',0)

                expected_numerical_integrand = (sp.sympify( str(self.cal_I) ) - sp.sympify('A + D*x1 + x0*B') ) * sp.sympify( str(self.exponentiated_monomial) )

                for summand in I_j_pole_part.summands:
                    np.testing.assert_array_equal(summand.factors[0].expolist, [(0,2,0,0)])
                    np.testing.assert_array_equal(summand.factors[0].coeffs, [1])
                    self.assertEqual(  sp.sympify( str(summand.factors[0].exponent) )  -  sp.sympify('-2 - eps0 - 3*eps1'),  0  )

            elif j == 1:
                expected_pole_part = sp.sympify('''
                                                     (  1/( 2*(-2 - eps0 - 3*eps1) + 0 + 1 ) * (A + B*x0 + C*x0**2*x1 + D*x1) +
                                                        1/( 2*(-2 - eps0 - 3*eps1) + 1 + 1 ) * (C*x0**2 + D) +
                                                        0  ) * x0**(-2 - eps0 - 3*eps1)
                                                ''').subs('x1',0)

                expected_numerical_integrand = (sp.sympify( str(self.cal_I) ) - sp.sympify('A + B*x0 + C*x0**2*x1 + D*x1')) * sp.sympify( str(self.exponentiated_monomial) )

            else:
                raise IndexError('`j` should only run over `(0,1)`!')

            self.assertEqual( (sp.sympify(str(I_j_pole_part)) - expected_pole_part).simplify() , 0)
            for summand in I_j_pole_part.summands:
                self.assertEqual(type(summand), Product)
                self.assertEqual(type(summand.factors[0]), ExponentiatedPolynomial)

            self.assertEqual(type(I_j_numerically_integrable_part), Product)
            self.assertEqual(type(I_j_numerically_integrable_part.factors[0]), ExponentiatedPolynomial)
            self.assertEqual( (sp.sympify(str(I_j_numerically_integrable_part)) - expected_numerical_integrand).simplify() , 0)

    def test_integrate_multiple_pole_parts(self):
        I_j_before = Product(self.exponentiated_monomial,self.regulator_poles,self.cal_I)
        I_j_after = PolynomialSum(*integrate_pole_part(I_j_before,0,1))
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
        I = Product(self.exponentiated_monomial,self.regulator_poles,self.cal_I)
        I_0 = integrate_pole_part(I,0)
        for i,term in enumerate(I_0):
            integrate_pole_part(term,1)
