from .expansion import *
from .expansion import _expand_singular_step, _expand_Taylor_step
from .algebra import Polynomial, ExponentiatedPolynomial, LogOfPolynomial, Product, Sum, _Expression
from .misc import flatten, sympify_expression
from nose.plugins.attrib import attr
import sympy as sp
import unittest

class TestSingularExpansion(unittest.TestCase):
    def setUp(self):
        self.unit_polynomial = Polynomial.from_expression('1', ['eps0','eps1'])
        self.p0 = Polynomial([(0,1),(1,0)], coeffs=[36, 12], polysymbols=['eps0','eps1'])
        self.p1 = Polynomial([(0,1),(1,0)], coeffs=[ 3,  1], polysymbols=['eps0','eps1'])
        self.p2 = Polynomial([(0,1)], coeffs=['1'], polysymbols=['eps0','eps1'])
        self.p3 = ExponentiatedPolynomial([(0,1),(1,0)], coeffs=[36, 12], polysymbols=['eps0','eps1'], exponent=-1)
        self.p4 = LogOfPolynomial([(0,1),(1,0)], coeffs=[36, 12], polysymbols=['eps0','eps1'])
        self.p5 = ExponentiatedPolynomial([(0,1),(1,0)], coeffs=[36, 12], polysymbols=['eps0','eps1'], exponent='eps0')
        self.p6_1 = Polynomial([(0,),(1,),(2,)], coeffs=['3/16','3/16','-9/16'], polysymbols=['eps'])
        self.p6_2 = ExponentiatedPolynomial([(2,),(3,),(4,)], coeffs=['9/16','-9/8','9/16'], polysymbols=['eps'], exponent=-1)
        self.p6 = Product(self.p6_1,self.p6_2)
        self.p7 = ExponentiatedPolynomial([(1,1),(1,2),(2,0),(2,1),(2,2),(3,0),(3,1),(4,0)],  
            coeffs=['3/4','-3/4','-3/16','-3/16','9/16','3/32','-9/32','9/256'],
            polysymbols=['n1','eps'], exponent=-1)

        self.numerator = self.unit_polynomial
        self.denominator = self.p0 * self.p1 * self.p2
        self.denominator = ExponentiatedPolynomial(self.denominator.expolist, self.denominator.coeffs, polysymbols=self.denominator.polysymbols, exponent=-1)
        self.rational_polynomial = Product(self.numerator, self.denominator)

    #@attr('active')
    def test_basic_checks(self):
        correct_input = Product(self.p0, self.p3)
        second_factor_wrong_exponent = Product(self.p0, self.p5)
        second_factor_wrong_type = Product(self.p0, self.p4)

        for expansion_function in (_expand_singular_step, expand_singular):
            # must have a rational polynomial (polynomial product of the form p * p**-1) in the first arg
            self.assertRaisesRegexp(TypeError, 'product.*must.*Product', expansion_function, self.p0, 0, 0)
            self.assertRaisesRegexp(TypeError, 'factor.*exponent.*-1', expansion_function, second_factor_wrong_exponent, 0, 0)
            self.assertRaisesRegexp(TypeError, '(A|a)ll.*factors.*Polynomial.*ExponentiatedPolynomial', expansion_function, second_factor_wrong_type, 0, 0)
            expansion_function(correct_input, 0, 0) # should not raise an error

        self.assertRaisesRegexp(AssertionError, 'indices.*orders.*same length', expand_singular, correct_input, indices=1, orders=[1,2])

    #@attr('active')
    def test_multiple_terms(self):
        unexpanded = Product(self.p3, self.p3, self.p2)

        expanded_0 = expand_singular(unexpanded, 0, 1)
        sympy_expanded_0 = sympify_expression(expanded_0)
        target_expanded_0 = sympify_expression('1/(1296*eps1) * eps0**0 - 1/(1944*eps1**2) * eps0**1')
        self.assertEqual( (sympy_expanded_0 - target_expanded_0).simplify(), 0 )

        expanded_0_1 = expand_singular(unexpanded, [0,1], [1,10])
        sympy_expanded_0_1 = sympify_expression(expanded_0_1)
        target_expanded_0_1 = target_expanded_0 # already a Laurent expansion
        self.assertEqual( (sympy_expanded_0_1 - target_expanded_0_1).simplify(), 0 )

    #@attr('active')
    def test_two_regulators_step(self):
        for error_t in (ValueError, OrderError):
            self.assertRaisesRegexp(error_t, 'lowest order.*higher than the requested order', _expand_singular_step, self.rational_polynomial, index=1, order=-2)

        expansion = _expand_singular_step(self.rational_polynomial, index=1, order=1)

        self.assertTrue(type(expansion) is Polynomial)
        self.assertEqual(expansion.number_of_variables, 2)
        for coeff in expansion.coeffs:
            self.assertTrue(type(coeff) is Product)

        # expansion in eps1 yields a simple pole --> expansion to order epsilon has three terms
        self.assertEqual(len(expansion.coeffs), 3)

        pole_order = sympify_expression('1/(12*eps0**2) * 1/eps1')
        constant_order = sympify_expression('-(12*3+36)*eps0/(12*12*eps0**4) * 1')
        order_epsilon = sympify_expression('9/(2*eps0**4) * eps1/2')

        target_expansion = pole_order + constant_order + order_epsilon
        self.assertEqual(target_expansion - sympify_expression(expansion), 0)

        # expand in the other regulator 'eps0'
        second_expansion = expansion.copy()
        for i, coeff in enumerate(expansion.coeffs):
            second_expansion.coeffs[i] = _expand_singular_step(coeff, index=0, order=0)

        # `target_expansion` is already expanded in 'eps0'
        self.assertEqual( (sympify_expression(expansion) - sympify_expression(second_expansion)).simplify() , 0)
        self.assertEqual( (target_expansion - sympify_expression(second_expansion)).simplify() , 0)

    #@attr('active')
    def test_high_level_function_one_regulator(self):
        for error_t in (ValueError, OrderError):
            self.assertRaisesRegexp(error_t, 'lowest order.*higher than the requested order', expand_singular, self.rational_polynomial, indices=1, orders=-2)

        expansion = expand_singular(self.rational_polynomial, indices=1, orders=1)

        self.assertTrue(type(expansion) is Polynomial)
        self.assertEqual(expansion.number_of_variables, 2)
        for coeff in expansion.coeffs:
            self.assertTrue(isinstance(coeff, _Expression))

        # expansion in eps1 yields a simple pole --> expansion to order epsilon has three terms
        self.assertEqual(len(expansion.coeffs), 3)

        pole_order = sympify_expression('1/(12*eps0**2) * 1/eps1')
        constant_order = sympify_expression('-(12*3+36)*eps0/(12*12*eps0**4) * 1')
        order_epsilon = sympify_expression('9/(2*eps0**4) * eps1/2')

        target_expansion = pole_order + constant_order + order_epsilon
        self.assertEqual(target_expansion - sympify_expression(expansion), 0)

    #@attr('active')
    def test_high_level_function_two_regulators(self):
        # expand in regulator 1 first, then in regulator 0
        expansion_1 = _expand_singular_step(self.rational_polynomial, index=1, order=1)
        expansion_1_0 = expansion_1.copy()
        for i, coeff in enumerate(expansion_1.coeffs):
            expansion_1_0.coeffs[i] = _expand_singular_step(coeff, index=0, order=0)
        flattened_expansion_1_0 = flatten(expansion_1_0, 1)

        # the high-level function should run exactly the commands above
        high_level_output = expand_singular(self.rational_polynomial, indices=[1,0], orders=[1,0])

        self.assertTrue(type(high_level_output) is Polynomial)
        for coeff in high_level_output.coeffs:
            self.assertTrue(isinstance(coeff, _Expression))

        self.assertEqual( (sympify_expression(high_level_output) - sympify_expression(flattened_expansion_1_0)).simplify() , 0)

    #@attr('active')
    def test_high_level_function_two_regulators_higher_order(self):
        expanded = expand_singular(self.rational_polynomial, indices=[1,0], orders=[3,2])
        target = sympify_expression('1/(12*eps0**2*eps1) - 1/(2*eps0**3) + 9*eps1/(4*eps0**4) -  9*eps1**2/eps0**5 + 135*eps1**3/(4*eps0**6)')
        self.assertEqual( (sympify_expression(expanded) - target).simplify() , 0)

    #@attr('active')
    def test_high_level_function_one_regulator_christoph_greub(self):
        expanded = expand_singular(self.p6, indices=[0], orders=[4])
        target = sympify_expression('1/(3*eps**2)+1/eps+2/3+eps/3-eps**3/3-(2*eps**4)/3')
        self.assertEqual( (sympify_expression(expanded) - target).simplify() , 0)

    #@attr('active')
    def test_high_level_function_two_regulators_christoph_greub(self):
        expanded = expand_singular(Product(self.p7, copy=False), indices=[0,1], orders=[0,4])
        target = sympify_expression('2/3+1/(3*eps**2)+1/eps+eps/3-eps**3/3-2*eps**4/3+4/(3*n1)+4/(3*eps*n1)+(4*eps)/(3*n1)+(4*eps**2)/(3*n1)+(4*eps**3)/(3*n1)+(4*eps**4)/(3*n1)')
        self.assertEqual( (sympify_expression(expanded) - target).simplify() , 0)

class TestTaylorExpansion(unittest.TestCase):
    def setUp(self):
        p0 = Polynomial.from_expression('5 + 2*x + 4*y + y**2', ['x','y'])
        self.p0 = ExponentiatedPolynomial(p0.expolist, p0.coeffs, polysymbols=p0.polysymbols, exponent=Polynomial.from_expression('2*x', ['x','y']))
        self.p1 = Polynomial.from_expression('3*x + y', ['x','y'])
        self.expression = Sum(self.p0, Product(self.p0, self.p1)) # (3*x + y)*(2*x + y**2 + 4*y)**(2*x) + (2*x + y**2 + 4*y)**(2*x)
        self.expected_expansion_in_x = sympify_expression('''
                                                         + x**0 / 0!  *  (y + 1)
                                                         + x**1 / 1!  *  (2*y*log(5 + y**2 + 4*y) + 2*log(5 + y**2 + 4*y) + 3) +
                                                         + x**2 / 2!  *  (y*(4*log(5 + y**2 + 4*y)**2 + 8/(5 + y**2 + 4*y)) + 4*log(5 + y**2 + 4*y)**2 + 12*log(5 + y**2 + 4*y) + 8/(5 + y**2 + 4*y))
                                                  ''')

    def test_error_messages(self):
        expression = Polynomial.from_expression('x', ['x','y'])
        self.assertRaisesRegexp(AssertionError, "order.*nonnegative integer", expand_Taylor, expression, indices=0, orders=-1)
        self.assertRaisesRegexp(AssertionError, "order.*nonnegative integer", expand_Taylor, expression, indices=0, orders=1.5)
        self.assertRaisesRegexp(IndexError, "out of bounds", expand_Taylor, expression, indices=4, orders=1)
        self.assertRaisesRegexp(AssertionError, 'indices.*orders.*same length', expand_Taylor, expression, indices=1, orders=[1,2])

    #@attr('active')
    def test_expand_Taylor_step(self):
        expansion_in_x = _expand_Taylor_step(self.expression, 0, 2)
        self.assertEqual( (sympify_expression(expansion_in_x) - self.expected_expansion_in_x).simplify() , 0)

    #@attr('active')
    def test_expand_Taylor(self):
        # the order of the expansion should not matter
        expansion_in_x_and_y = expand_Taylor(self.expression, [0,1], [2,0])
        expansion_in_y_and_x = expand_Taylor(self.expression, [1,0], [0,2])
        expected_expansion = self.expected_expansion_in_x.subs('y',0) # zeroth order in y

        self.assertEqual( (sympify_expression(expansion_in_x_and_y) - expected_expansion).simplify() , 0)
        self.assertEqual( (sympify_expression(expansion_in_y_and_x) - expected_expansion).simplify() , 0)

#@attr('active')
class TestExpandSympy(unittest.TestCase):
    #@attr('active')
    def test_error_messages(self):
        expression = 'a + b'
        variables = ['a', 'b']
        orders = [1,3]

        self.assertRaisesRegexp(AssertionError, '(N|n)umber of variables \(2\).*must.*equal.*number of orders \(3\)', expand_sympy, expression, variables, [0,0,1])
        self.assertRaisesRegexp(AssertionError, 'variables.*must.*symbols', expand_sympy, expression, ['a+b','x'], orders)
        self.assertRaisesRegexp(AssertionError, 'orders.*must.*vector', expand_sympy, expression, variables, [[0,1],[1,2]])

        expand_sympy(expression, variables, orders) # should be ok

    #@attr('active')
    def test_requested_order_too_low(self):
        expression = 'x**3'
        variables = ['x']
        orders = [1]
        for error_t in (ValueError, OrderError):
            self.assertRaisesRegexp(error_t, 'lowest order.*x.*\(3\).*higher than.*requested.*\(1\)', expand_sympy, expression, variables, orders)

    #@attr('active')
    def test_1d(self):
        expression = 'exp(x)/x'
        variables = ['x']
        orders = [1]

        poly = expand_sympy(expression, variables, orders)

        self.assertTrue( type(poly) is Polynomial )
        np.testing.assert_array_equal(poly.expolist,            [[-1],[0],[  1  ]])
        np.testing.assert_array_equal(poly.coeffs,    sympify_expression([ 1 , 1 , '1/2']))

    #@attr('active')
    def test_2d(self):
        expression = '1/(eps+alpha)'
        variables = sympify_expression(['alpha', 'eps'])
        orders = [0,1]

        # expansion in 'alpha' first
        alpha_first = expand_sympy(expression, variables, orders)
        target_alpha_first = sympify_expression('1/eps')

        self.assertEqual(  (sympify_expression(alpha_first) - target_alpha_first).simplify() , 0  )

        self.assertTrue( type(alpha_first) is Polynomial )
        self.assertEqual(alpha_first.polysymbols, variables)
        np.testing.assert_array_equal(alpha_first.expolist, [[0,0]])

        alpha_poly = alpha_first.coeffs[0]
        self.assertTrue( type(alpha_poly) is Polynomial )
        np.testing.assert_array_equal( alpha_poly.expolist, [[0,-1]] )
        np.testing.assert_array_equal( alpha_poly.coeffs, [1] )

        # expansion in 'eps' first
        variables = sympify_expression(['eps', 'alpha'])
        orders = [1,0]
        eps_first = expand_sympy(expression, variables, orders)
        target_eps_first = sympify_expression('1/alpha - eps/alpha**2')

        self.assertEqual(  (sympify_expression(eps_first) - target_eps_first).simplify() , 0  )

        self.assertTrue( type(eps_first) is Polynomial )
        self.assertEqual(eps_first.polysymbols, variables)
        np.testing.assert_array_equal(eps_first.expolist, [[0,0],[1,0]])

        eps_to_the_0 = eps_first.coeffs[0]
        eps_to_the_1 = eps_first.coeffs[1]
        self.assertTrue( type(eps_to_the_0) is Polynomial )
        self.assertTrue( type(eps_to_the_1) is Polynomial )
        np.testing.assert_array_equal( eps_to_the_0.expolist, [[0,-1]] )
        np.testing.assert_array_equal( eps_to_the_1.expolist, [[0,-2]] )
        np.testing.assert_array_equal( eps_to_the_0.coeffs, [ 1] )
        np.testing.assert_array_equal( eps_to_the_1.coeffs, [-1] )

    #@attr('active')
    def test_truncation_field(self):
        expression = 'x + x**2'
        variables = ['x']

        expansion = expand_sympy(expression, variables, orders=[1])
        target_expansion = Polynomial([[1]],[1],['x'])
        self.assertEqual(  (sympify_expression(expansion) - sympify_expression(target_expansion)).simplify() , 0  )
        self.assertTrue(expansion.truncated is True)

        expansion = expand_sympy(expression, variables, orders=[4])
        target_expansion = Polynomial([[1],[2]],[1,1],['x'])
        self.assertEqual(  (sympify_expression(expansion) - sympify_expression(target_expansion)).simplify() , 0  )
        self.assertTrue(expansion.truncated is False)

    #@attr('active')
    def test_missing_intermediate_order(self):
        expression = 'exp(3*EulerGamma*eps)*gamma(3*eps)'
        variables = ['eps']

        expansion = expand_sympy(expression, variables, orders=[4])

        target_expansion_expolist = np.arange(6).reshape([6,1]) - 1
        target_expansion_coeffs = [
                                      '1/(3)',                                         # eps ** -1
                                      0,                                               # eps **  0
                                      'pi**2/4',                                       # eps **  1
                                      '3*polygamma(2,1)/2',                            # eps **  2
                                      '27*pi**4/160',                                  # eps **  3
                                      '9/40*(5*pi**2*polygamma(2,1)+3*polygamma(4,1))' # eps **  4
                                  ]
        target_expansion = Polynomial(target_expansion_expolist, target_expansion_coeffs, ['eps'])

        self.assertEqual( sympify_expression(expansion - target_expansion).simplify() , 0 )

        np.testing.assert_array_equal(expansion.expolist, target_expansion.expolist)
        np.testing.assert_array_equal(expansion.coeffs, target_expansion.coeffs)
        self.assertTrue(expansion.truncated is True)

        for coeff in expansion.coeffs:
            self.assertTrue(isinstance(coeff, sp.Expr))

    #@attr('active')
    def test_higher_pole(self):
        expression = 'gamma(eps)/eps'
        variables = ['eps']

        expansion = expand_sympy(expression, variables, orders=[0])

        target_expansion_expolist = np.arange(3).reshape([3,1]) - 2
        target_expansion_coeffs = [
                                      '1',                         # eps ** -2
                                      '-EulerGamma',               # eps ** -1
                                      'EulerGamma**2/2 + pi**2/12' # eps **  0
                                  ]
        target_expansion = Polynomial(target_expansion_expolist, target_expansion_coeffs, ['eps'])

        self.assertEqual( sympify_expression(expansion - target_expansion).simplify() , 0 )

        np.testing.assert_array_equal(expansion.expolist, target_expansion.expolist)
        np.testing.assert_array_equal(expansion.coeffs, target_expansion.coeffs)
        self.assertTrue(expansion.truncated is True)

        for coeff in expansion.coeffs:
            self.assertTrue(isinstance(coeff, sp.Expr))

    #@attr('active')
    def test_nontrivial_higher_pole(self):
        expression = 'gamma(eps+2)/eps^2 + a/eps^2'
        variables = ['eps']

        expansion = expand_sympy(expression, variables, orders=[1])

        target_expansion_expolist = np.arange(4).reshape([4,1]) - 2
        target_expansion_coeffs = [
                                      '1 + a',                                                                                        # eps ** -2
                                      '-EulerGamma+1',                                                                                # eps ** -1
                                      'pi**2/12 + EulerGamma**2/2 - EulerGamma',                                                      # eps **  0
                                      '-EulerGamma*pi**2/12 - 1/3 + polygamma(2, 2)/6 - EulerGamma**3/6 + EulerGamma**2/2 + pi**2/12' # eps **  1
                                  ]
        target_expansion = Polynomial(target_expansion_expolist, target_expansion_coeffs, ['eps'])

        self.assertEqual( sympify_expression(expansion - target_expansion).simplify() , 0 )

        np.testing.assert_array_equal(expansion.expolist, target_expansion.expolist)
        np.testing.assert_array_equal(expansion.coeffs, target_expansion.coeffs)
        self.assertTrue(expansion.truncated is True)

        for coeff in expansion.coeffs:
            self.assertTrue(isinstance(coeff, sp.Expr))
