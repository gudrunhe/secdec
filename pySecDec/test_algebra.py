from .algebra import *
from .algebra import _Expression
from .misc import sympify_expression
import sympy as sp
import unittest
from nose.plugins.attrib import attr

#@attr('active')
class TestFunction(unittest.TestCase):
    def setUp(self):
        self.polysymbols = ['x0', 'x1', 'x2', 'x3']
        self.arg0 = Polynomial.from_expression('x0', self.polysymbols)
        self.arg1 = Polynomial.from_expression('x1', self.polysymbols)
        self.simplifyable_arg = Polynomial([[0,0,1,1],[0,0,1,1],[1,0,0,0]], [1,1,1], self.polysymbols)
        self.f = Function('f', self.arg0, self.arg1)
        self.g = Function('g', self.arg0, self.arg1, self.simplifyable_arg)

    #@attr('active')
    def test_init(self):
        # do not use variables defined in `setUp` here because of the last test in this function
        polysymbols = ['x0', 'x1', 'x2', 'x3']
        arg0 = Polynomial.from_expression('x0', polysymbols)
        arg1 = Polynomial.from_expression('x1', polysymbols)
        wrong_arg = Polynomial.from_expression('x4', polysymbols + ['x4'])
        self.assertRaisesRegexp(AssertionError, 'same.*number of variables.*all arguments', Function, 'f', arg0, arg1, wrong_arg)
        f = Function('f', arg0, arg1)

        self.assertEqual(f.number_of_variables, len(polysymbols))
        self.assertEqual( str(f) , 'f( + (1)*x0, + (1)*x1)' )

        # made a copy?
        arg0.expolist += 1
        self.assertEqual( str(f) , 'f( + (1)*x0, + (1)*x1)' )

    #@attr('active')
    def test_copy(self):
        f = self.f.copy()
        f.arguments[0].expolist[:,0] = 0
        f.arguments[0].expolist[:,1] = 1
        self.assertNotEqual( str(self.f) , str(f) )
        self.assertEqual( str(f) , 'f( + (1)*x1, + (1)*x1)' )

    #@attr('active')
    def test_simplify(self):
        unsimplified_g = self.g.copy()
        str_unsimplified_g = str(unsimplified_g)

        simplified_g = unsimplified_g.simplify()
        str_simplified_g = str(simplified_g)

        self.assertNotEqual(str_simplified_g, str_unsimplified_g)
        self.assertEqual( (sympify_expression(str_simplified_g) - sympify_expression(str_unsimplified_g)).simplify() , 0 )
        self.assertEqual( str(self.simplifyable_arg.simplify()) , str(simplified_g.arguments[-1]) )

    #@attr('active')
    def test_derive(self):
        # simple derivative
        dfd0 = self.f.derive(0)
        self.assertEqual( (sympify_expression(dfd0) - sympify_expression('dfd0(x0,x1)')).simplify() , 0 )

        # g( + (1)*x0, + (1)*x1, + (1)*x2*x3 + (1)*x2*x3 + (1)*x0)

        # derivatives with chain rule
        dgd2 = self.g.derive(2)
        self.assertEqual( (sympify_expression(dgd2) - sympify_expression('dgd2(x0,x1,2*x2*x3+x0)*2*x3')).simplify() , 0 )

        dgd0 = self.g.derive(0)
        self.assertEqual( (sympify_expression(dgd0) - sympify_expression('dgd0(x0,x1,2*x2*x3+x0) + dgd2(x0,x1,2*x2*x3+x0)')).simplify() , 0 )

    #@attr('active')
    def test_derivative_symbols(self):
        polysymbols = ['x','y','z']
        x = Polynomial.from_expression('x', polysymbols)
        y = Polynomial.from_expression('y', polysymbols)
        z = Polynomial.from_expression('z', polysymbols)

        derivatives = set()
        f = Function('f', x, y, derivative_symbols=derivatives).replace(-1,0)

        self.assertEqual(derivatives, set(['f']))

        f.derive(0)
        self.assertEqual(derivatives, set(['f','dfd0']))

        f.derive(2)
        self.assertEqual(derivatives, set(['f','dfd0']))

        f.derive(1).derive(0)
        self.assertEqual(derivatives, set(['f','dfd0','dfd1','ddfd0d1']))

    #@attr('active')
    def test_derivative_sorting(self):
        polysymbols = ['x','y','z']
        x = Polynomial.from_expression('x', polysymbols)
        y = Polynomial.from_expression('y', polysymbols)
        z = Polynomial.from_expression('z', polysymbols)

        derivatives = set()
        f = Function('f', x, y, z, derivative_symbols=derivatives)

        target_derivatives = set(['f'])
        self.assertEqual(derivatives, target_derivatives)

        dfd0d1 = f.derive(1).derive(0).simplify()
        target_derivatives.update(['f','dfd1','ddfd0d1'])
        self.assertEqual(  (sympify_expression(dfd0d1) - sympify_expression('ddfd0d1(x,y,z)')).simplify()  ,  0  )
        self.assertEqual(derivatives, target_derivatives)

        dddfd0d0d1 = f.derive(1).derive(0).derive(0).simplify()
        target_derivatives.update(['f','dfd1','ddfd0d1','dddfd0d0d1'])
        self.assertEqual(  (sympify_expression(dddfd0d0d1) - sympify_expression('dddfd0d0d1(x,y,z)')).simplify()  ,  0  )
        self.assertEqual(derivatives, target_derivatives)

        ddddfd0d0d1d2_with_copy = f.derive(1).copy().derive(0).simplify().derive(2).copy().derive(0).simplify()
        target_derivatives.update(['f','dfd1','ddfd0d1','dddfd0d1d2','ddddfd0d0d1d2'])
        self.assertEqual(  (sympify_expression(ddddfd0d0d1d2_with_copy) - sympify_expression('ddddfd0d0d1d2(x,y,z)')).simplify()  ,  0  )
        self.assertEqual(derivatives, target_derivatives)

    #@attr('active')
    def test_derivative_sorting_composite_arg(self):
        polysymbols = ['x','y','z']
        x,y,z = (Polynomial.from_expression(symbol, polysymbols) for symbol in polysymbols)

        derivatives = set()
        f = Function('f', x-y, x+y, z, derivative_symbols=derivatives)

        target_derivatives = set(['f'])
        self.assertEqual(derivatives, target_derivatives)

        ddfd0d1 = f.derive(1).derive(0).simplify()
        target_derivatives.update(['f','dfd1','dfd0','ddfd0d1','ddfd0d0','ddfd1d1'])
        self.assertEqual(  (sympify_expression(ddfd0d1) - sympify_expression('ddfd1d1(x-y,x+y,z)-ddfd0d0(x-y,x+y,z)')).simplify()  ,  0  )
        self.assertEqual(derivatives, target_derivatives)

    #@attr('active')
    def test_derivative_tracking(self):
        x = Polynomial.from_expression('x', ['x','y'])
        y = Polynomial.from_expression('y', ['x','y'])
        poly = x**2*y + y**2
        func = Function('f', x, y)
        func.derive(0).derive(1)
        target_derivative_tracks = {(1,0): [0,(0,0)], (1,1) : [1,(1,0)]}
        self.assertEqual(func.derivative_tracks, target_derivative_tracks)

    #@attr('active')
    def test_derivative_tracking_with_replace(self):
        x = Polynomial.from_expression('x', ['x','y','z'])
        y = Polynomial.from_expression('y', ['x','y','z'])
        z = Polynomial.from_expression('z', ['x','y','z'])
        poly = Polynomial.from_expression('x**2*y + y**2 + z', ['x','y','z'])
        func = Function('f', x, y, z)
        func.derive(0).replace(0,0).derive(1).replace(1,0)
        target_derivative_tracks = {(1,0,0): [0,(0,0,0)], (1,1,0) : [1,(1,0,0)]}
        self.assertEqual(func.derivative_tracks, target_derivative_tracks)

    #@attr('active')
    def test_str(self):
        # Note: ``repr`` is checked by the doctests
        x, y, z = (Polynomial.from_expression(symbol, ['x','y','z']) for symbol in ['x','y','z'])
        func = Function('func', x, y, z)
        self.assertEqual(str(func), 'func( + (1)*x, + (1)*y, + (1)*z)')

    #@attr('active')
    def test_compute_own_derivatives(self):
        polysymbols = ['x','y']
        x,y = (Polynomial.from_expression(symbol, polysymbols) for symbol in polysymbols)

        func = Function('F', x**2*y, y**2)
        func.derive(0).derive(1)
        func.derive(1).derive(0).derive(0)

        target_derivatives_as_strings = \
        {
            (1,0) : 'dFd0(x**2*y,y**2) * 2*x*y',
            (0,1) : 'dFd0(x**2*y,y**2) * x**2 + dFd1(x**2*y,y**2) * 2*y',
            (1,1) : 'dFd0(x**2*y,y**2) * 2*x + ddFd0d0(x**2*y,y**2) * x**2*2*x*y + ddFd0d1(x**2*y,y**2) * 2*y*2*x*y',
            (2,1) : 'dFd0(x**2*y,y**2) * 2 + ddFd0d0(x**2*y,y**2) * 2*x*2*x*y + ddFd0d0(x**2*y,y**2) * 6*x**2*y + dddFd0d0d0(x**2*y,y**2) * x**2*2*x*y*2*x*y +' + \
                    'ddFd0d1(x**2*y,y**2) * 2*y*2*y + dddFd0d0d1(x**2*y,y**2) * 2*y*2*x*y*2*x*y',
            (2,0) : 'dFd0(x**2*y,y**2) * 2*y + ddFd0d0(x**2*y,y**2) * (2*x*y)**2',
            (3,0) : 'ddFd0d0(x**2*y,y**2) * 2*y*2*x*y + ddFd0d0(x**2*y,y**2) * 2*4*x**1*y**2 + dddFd0d0d0(x**2*y,y**2) * 4*x**2*y**2*2*x*y'
        }
        target_derivatives = {key : sympify_expression(value).expand()
                              for key, value in target_derivatives_as_strings.items()}

        recomputed_derivatives_as_expressions = func.compute_derivatives()
        recomputed_derivatives = {key : sympify_expression(value).expand()
                                  for key, value in recomputed_derivatives_as_expressions.items()}

        self.assertEqual(recomputed_derivatives, target_derivatives)

    #@attr('active')
    def test_compute_other_derivatives(self):
        polysymbols = ['x','y']
        x,y = (Polynomial.from_expression(symbol, polysymbols) for symbol in polysymbols)
        poly = Polynomial.from_expression('x**2*y + y**2', polysymbols)

        func = Function('f', x, y)
        func.derive(0).derive(1)
        func.derive(1).derive(0).derive(0)

        other_expression = x*y

        target_derivatives_as_strings = \
        {
            (1,0) : 'y',
            (0,1) : 'x',
            (1,1) : '1',
            (2,1) : '0'
        }
        target_derivatives = {key : sympify_expression(value)
                              for key, value in target_derivatives_as_strings.items()}

        recomputed_derivatives_as_expressions = func.compute_derivatives(other_expression)
        recomputed_derivatives = {key : sympify_expression(value)
                                  for key, value in recomputed_derivatives_as_expressions.items()}

        self.assertEqual(recomputed_derivatives, target_derivatives)

class TestPolynomial(unittest.TestCase):
    def test_init(self):
        # Proper instantiation
        Polynomial([(0,1),(1,0),(2,1)],['A','B','C'])

        # exponents should be integer
        self.assertRaisesRegexp(TypeError, "(e|E)xpolist.*integer", Polynomial, [(0,1.5),(1,0),(2,1)], ['A','B','C'])

        # Length mismatch between coeffs and expolist
        self.assertRaisesRegexp(AssertionError, "same length", Polynomial, [(0,1),(1,0),(2,1)], ['A','B','C','D'])

        # entries of expolist have variable length
        self.assertRaisesRegexp(AssertionError, "expolist.*same length", Polynomial, [(0,1,2),(1,0),(2,1)], ['A','B','C'])

        # same number of variables in coeff if coeffs are `_Expression`s
        self.assertRaisesRegexp(AssertionError, "same number of variables.*(p|P)olynomial.*coeffs", Polynomial, [(0,2),(1,0),(2,1)], [Polynomial.from_expression('1', ['z']),'B','C'])

        Polynomial([(0,2),(1,0),(2,1)], [Polynomial.from_expression('1', ['x0','x1']),'B','C']) # OK

    def test_string_form(self):
        polynomial1 = Polynomial([(0,1),(1,0),(2,1),(0,0)],['A','B','C','D'])
        string_polynomial1 = " + (A)*x1 + (B)*x0 + (C)*x0**2*x1 + (D)"
        self.assertEqual(str(polynomial1), string_polynomial1)
        self.assertEqual(repr(polynomial1), string_polynomial1)

        polynomial2 = Polynomial([(0,1),(1,0),(2,1),(0,0)],['A','B','C','D'], polysymbols='z')
        string_polynomial2 = string_polynomial1.replace('x','z')
        self.assertEqual(str(polynomial2), string_polynomial2)
        self.assertEqual(repr(polynomial2), string_polynomial2)

        polynomial3 = Polynomial([(0,1),(1,0),(2,1),(0,0)],['A','B','C','D'], polysymbols=['x','y'])
        string_polynomial3 = string_polynomial1.replace('x0','x').replace('x1','y')
        self.assertEqual(str(polynomial3), string_polynomial3)
        self.assertEqual(repr(polynomial3), string_polynomial3)

    def test_copy(self):
        A, Afoo = sp.symbols('A Afoo')

        polynomial1 = Polynomial([(0,1),(1,0),(2,1),(0,0)],['A','B','C','D'])
        polynomial2 = polynomial1.copy()

        self.assertEqual(str(polynomial1), str(polynomial2))

        polynomial1.expolist[0,0] = 5
        self.assertEqual(polynomial1.expolist[0,0],5)
        self.assertEqual(polynomial2.expolist[0,0],0)

        polynomial1.coeffs[0] = Afoo
        self.assertEqual(polynomial1.coeffs[0],Afoo)
        self.assertEqual(polynomial2.coeffs[0],A)

    #@attr('active')
    def test_has_constant_term(self):
        self.assertTrue(Polynomial([(0,1),(1,0),(2,1),(0,0)],['A','B','C','D']).has_constant_term())
        self.assertTrue(Polynomial([(0,1),(0,0),(1,0),(2,1),(4,0)],['A','B','C','D',1]).has_constant_term())

        self.assertFalse(Polynomial([(0,1),(1,0),(2,1),(4,0)],['A','B','C','D']).has_constant_term())
        self.assertFalse(Polynomial([(0,1),(2,1),(4,0)],['A','B','D']).has_constant_term())

        self.assertTrue(Polynomial([(0,1,1),(1,0,2),(2,1,3),(0,0,4)],['A','B','C','D']).has_constant_term([0,1]))
        self.assertTrue(Polynomial([(0,2,1),(1,3,0),(2,4,1),(0,5,0)],['A','B','C','D']).has_constant_term([0,2]))
        self.assertFalse(Polynomial([(0,1,1),(1,0,2),(2,1,3),(0,0,4)],['A','B','C','D']).has_constant_term())

    def test_becomes_zero_for(self):
        self.assertTrue(Polynomial([(0,1,1,0),(2,1,0,5)],['A','B']).becomes_zero_for([1,0]))
        self.assertTrue(Polynomial([(0,1),(0,1),(1,5),(0,1),(4,4)],['A','B','C','D',1]).becomes_zero_for([1]))

        self.assertFalse(Polynomial([(0,1),(1,0),(2,1),(4,0)],['A','B','C','D']).becomes_zero_for([1]))
        self.assertFalse(Polynomial([(0,1,0,1),(2,1,0,5),(0,0,3,5)],['A','B','C']).becomes_zero_for([1,0]))

    def test_derive(self):
        polynomial = Polynomial([(2,1),(0,1)],['A', 'B'])

        derivative_0 = sympify_expression( str(polynomial.derive(0)) )
        target_derivative_0 = sympify_expression('2*A*x0*x1')
        self.assertEqual( (derivative_0 - target_derivative_0).simplify() , 0 )

        derivative_1 = sympify_expression( str(polynomial.derive(1)) )
        target_derivative_1 = sympify_expression('A*x0**2 + B')
        self.assertEqual( (derivative_1 - target_derivative_1).simplify() , 0 )

    #@attr('active')
    def test_derive_with_coeffs(self):
        # A*x0 * x0**2*x1 + B * x1 = A * x0**3*x1 + B * x1
        expr = Polynomial([(2,1),(0,1)],[Polynomial([(1,0)], ['A']), 'B'])

        derivative_0 = sympify_expression(expr.derive(0))
        target_derivative_0 = sympify_expression('3*A * x0**2*x1')

        self.assertEqual( (derivative_0 - target_derivative_0).simplify() , 0)

    #@attr('active')
    def test_simplify(self):
        expr1 = Polynomial([(2,1),(0,1)],[Polynomial([(0,0)], [0]), 'B']).simplify()
        expr2 = Polynomial([(2,1),(0,1)],['A', Polynomial([(5,2)], [0])]).simplify()

        np.testing.assert_array_equal(expr1.expolist, [[0,1]])
        np.testing.assert_array_equal(expr1.coeffs, [sp.symbols('B')])

        np.testing.assert_array_equal(expr2.expolist, [[2,1]])
        np.testing.assert_array_equal(expr2.coeffs, [sp.symbols('A')])

    #@attr('active')
    def test_refactorize(self):
        prod = Product(Polynomial([(0,0,0)],[1],'t'), Polynomial([(1,1,0),(1,0,1)],["-s12","-s23"],'t'))

        self.assertEqual(str(prod.factors[0]), ' + (1)')
        self.assertEqual(str(prod.factors[1]), ' + (-s12)*t0*t1 + (-s23)*t0*t2')

        copy0 = prod.copy()
        copy1 = prod.copy()
        copy2 = prod.copy()

        refactorize(copy0) # refactorize all parameters -> should find the factorization of parameter 0
        refactorize(copy1,0) # refactorize parameter 0 -> should find a factorization
        refactorize(copy2,1) # refactorize parameter 1 -> should NOT find the factorization of parameter 0

        self.assertEqual(str(copy0.factors[0]), str(copy1.factors[0]))
        self.assertEqual(str(copy0.factors[1]), str(copy1.factors[1]))

        self.assertEqual(str(copy1.factors[0]), ' + (1)*t0')
        self.assertEqual(str(copy1.factors[1]), ' + (-s12)*t1 + (-s23)*t2')

        self.assertEqual(str(copy2.factors[0]), ' + (1)')
        self.assertEqual(str(copy2.factors[1]), ' + (-s12)*t0*t1 + (-s23)*t0*t2')

    #@attr('active')
    def test_refactorize_multiple_parameters(self):
        prod = Product(Polynomial([(0,0,0)],[1],'t'), Polynomial([(1,1,3),(1,0,1)],["-s12","-s23"],'t'))

        self.assertEqual(str(prod.factors[0]), ' + (1)')
        self.assertEqual(str(prod.factors[1]), ' + (-s12)*t0*t1*t2**3 + (-s23)*t0*t2')

        copy0 = prod.copy()
        copy1 = prod.copy()

        refactorize(copy0,0,2) # refactorize parameter 0 and 2 -> should find a factorization for both
        refactorize(copy1,0,1) # refactorize parameter 0 and 1 -> should only find the factorization of parameter 0

        self.assertEqual(str(copy0.factors[0]), ' + (1)*t0*t2')
        self.assertEqual(str(copy0.factors[1]), ' + (-s12)*t1*t2**2 + (-s23)')

        self.assertEqual(str(copy1.factors[0]), ' + (1)*t0')
        self.assertEqual(str(copy1.factors[1]), ' + (-s12)*t1*t2**3 + (-s23)*t2')

    #@attr('active')
    def test_refactorize_polynomial(self):
        poly = Polynomial([(1,1,0),(1,0,1)],["-s12","-s23"],'t')

        copy0 = poly.copy()
        copy1 = poly.copy()
        copy2 = poly.copy()

        copy0 = copy0.refactorize() # refactorize all parameters -> should find the factorization of parameter 0
        copy1 = copy1.refactorize(0) # refactorize parameter 0 -> should find a factorization
        copy2 = copy2.refactorize(1) # refactorize parameter 1 -> should NOT find the factorization of parameter 0

        self.assertEqual(str(copy0.factors[0]), str(copy1.factors[0]))
        self.assertEqual(str(copy0.factors[1]), str(copy1.factors[1]))

        self.assertEqual(str(copy1.factors[0]), ' + (1)*t0')
        self.assertEqual(str(copy1.factors[1]), ' + (-s12)*t1 + (-s23)*t2')

        self.assertEqual(str(copy2.factors[0]), ' + (1)')
        self.assertEqual(str(copy2.factors[1]), ' + (-s12)*t0*t1 + (-s23)*t0*t2')

    #@attr('active')
    def test_refactorize_exponentiated_polynomial(self):
        exponentiated_poly = ExponentiatedPolynomial([(1,1,0),(1,0,1)],["-s12","-s23"],3,'t')

        copy0 = exponentiated_poly.copy()
        copy1 = exponentiated_poly.copy()
        copy2 = exponentiated_poly.copy()

        copy0 = copy0.refactorize() # refactorize all parameters -> should find the factorization of parameter 0
        copy1 = copy1.refactorize(0) # refactorize parameter 0 -> should find a factorization
        copy2 = copy2.refactorize(1) # refactorize parameter 1 -> should NOT find the factorization of parameter 0

        self.assertEqual(sympify_expression(copy0.factors[0] - copy1.factors[0]).simplify(), 0)
        self.assertEqual(sympify_expression(copy0.factors[1] - copy1.factors[1]).simplify(), 0)

        self.assertEqual(sympify_expression(str(copy1.factors[0])+ ' - ' + ' (+(1)*t0)**(3)').simplify(), 0)
        self.assertEqual(sympify_expression(str(copy1.factors[1])+ ' - ' + ' + ((-s12)*t1 + (-s23)*t2)**3').simplify(), 0)

        self.assertEqual(sympify_expression(str(copy2.factors[0])+ ' - ' + ' 1').simplify(), 0)
        self.assertEqual(sympify_expression(str(copy2.factors[1])+ ' - ' + ' ((-s12)*t0*t1 + (-s23)*t0*t2)**3').simplify(), 0)

    #@attr('active')
    def test_simplify_to_zero(self):
        zero = Polynomial([[0,0,0]]*11,[-1,0,1,0,0,-2,0,0,0,1,1]).simplify()

        np.testing.assert_array_equal(zero.coeffs, [0])
        np.testing.assert_array_equal(zero.expolist, [[0,0,0]])

class TestFormerSNCPolynomial(unittest.TestCase):
    def setUp(self):
        self.p0 = Polynomial([(0,1),(1,0),(1,1)],[1,1,3])
        self.p1 = Polynomial([(1,1),(1,1),(2,1)],[1,1,2])
        self.p2 = Polynomial([(0,0),(1,0),(1,1)],[5,6,7])

    def test_init(self):
        self.assertRaisesRegexp(AssertionError, "coeffs.*one.dimensional", Polynomial, [(0,3),(1,2)],[[1,4],[2,3]])
        Polynomial([(1,2),(1,0),(2,1)], [1,2,3])

    def test_simplify(self):
        polynomial = Polynomial([[1,2],[1,2],[2,2]  ,  [2,1],[2,1],[3,1]  ,  [2,2],[2,2],[3,2]], [1,1,2  ,  1,1,2  ,  3,3,6])
        polynomial.simplify()
        self.assertEqual( (sympify_expression(str(polynomial)) - sympify_expression(" + (2)*x0*x1**2 + (8)*x0**2*x1**2 + (2)*x0**2*x1 + (2)*x0**3*x1 + (6)*x0**3*x1**2")).simplify() , 0)

        # should have minimal number of terms
        self.assertEqual(len(polynomial.coeffs), 5)
        self.assertEqual(polynomial.expolist.shape, (5,2))

    def test_mul(self):
        self.assertRaisesRegexp(AssertionError, "Number of varibales must be equal for both factors in \*", lambda: self.p0 * Polynomial([(0,0,3)],[4]))

        poly = Polynomial([(0,2),(1,0)],[1,'C'])
        intprod = 5 * poly
        self.assertEqual( (sympify_expression(str(intprod)) - sympify_expression(' + (5)*x1**2 + (5*C)*x0')).simplify() , 0)
        # should have minimal number of terms
        self.assertEqual(len(intprod.coeffs), 2)
        self.assertEqual(intprod.expolist.shape, (2,2))

        prod = self.p0 * self.p1

        #                              expolist    coefficient
        target_prod_expolist = np.array([[1,2],   #     1
        #                                [1,2],   #     1
                                         [2,2],   #     2
                                         [2,1],   #     1
        #                                [2,1],   #     1
                                         [3,1],   #     2
        #                                [2,2],   #     3
        #                                [2,2],   #     3
                                         [3,2]])  #     6

        target_prod_coeffs = np.array([2,  8  ,  2,  2  ,      6])


        # the order of the terms does not matter but for array comparison, it must be fixed
        sortkey_target = argsort_2D_array(target_prod_expolist)
        sorted_target_prod_expolist = target_prod_expolist[sortkey_target]
        sorted_target_prod_coeffs = target_prod_coeffs[sortkey_target]

        sortkey_prod = argsort_2D_array(prod.expolist)
        sorted_prod_expolist = prod.expolist[sortkey_prod]
        sorted_prod_coeffs = prod.coeffs[sortkey_prod]

        np.testing.assert_array_equal(sorted_prod_expolist, sorted_target_prod_expolist)
        np.testing.assert_array_equal(sorted_prod_coeffs, sorted_target_prod_coeffs)

    def test_add(self):
        self.assertRaisesRegexp(AssertionError, "Number of varibales must be equal for both polynomials in \+", lambda: self.p0 + Polynomial([(0,0,3)],[4]))

        self.assertEqual( str(10 + Polynomial([(1,0)],[1])) , ' + (10) + (1)*x0' )
        self.assertEqual( str(10 + Polynomial([(0,0),(1,0)],[2,1])) , ' + (12) + (1)*x0' )

        polysum = self.p0 + self.p2
        self.assertEqual( (sympify_expression(str(polysum)) - sympify_expression(" + (1)*x1 + (7)*x0 + (10)*x0*x1 + (5)")).simplify() , 0)

        # should have minimal number of terms
        self.assertEqual(len(polysum.coeffs), 4)
        self.assertEqual(polysum.expolist.shape, (4,2))


    def test_negation(self):
        neg_p2 = - self.p2
        np.testing.assert_array_equal(neg_p2.coeffs, -self.p2.coeffs)
        np.testing.assert_array_equal(neg_p2.expolist, self.p2.expolist)

    def test_subtract(self):
        polysum = self.p0 + self.p2
        polysum -= self.p2
        self.assertEqual(str(polysum),str(self.p0))

    #@attr('active')
    def test_pow_sympy(self):
        poly = Polynomial([[1]], [1])
        poly_pow_poly = poly**poly
        expected_string_form_1 = '( + (1)*x0) ** ( + (1)*x0)'
        self.assertEqual(str(poly_pow_poly), expected_string_form_1)

        poly_pow_string = poly ** 'a+5*b'
        self.assertTrue(isinstance(poly_pow_string.exponent.coeffs[0], sp.Expr)) # ``a+5*b`` should be converted to sympy expression
        poly_pow_sympy = poly ** sympify_expression('a+5*b')
        expected_string_form_2 = '( + (1)*x0) ** ( + (a + 5*b))'
        for i, poly_pow in enumerate((poly_pow_string, poly_pow_sympy)):
            print(i)
            self.assertEqual(str(poly_pow), expected_string_form_2)
            self.assertTrue(type(poly_pow) is Pow)

        poly_pow_negative_int = poly ** -1
        expected_string_form_3 = '( + (1)*x0) ** ( + (-1))'
        self.assertEqual(str(poly_pow_negative_int), expected_string_form_3)
        self.assertTrue(isinstance(poly_pow_negative_int.exponent.coeffs[0], sp.Expr)) # ``-1`` should be converted to sympy expression

        poly_pow_float = poly ** -1.5
        # expect something like ``( + (1)*x0) ** ( + (-1.50000000000000))``
        # the exact representation of the floating point number can be system dependent --> use a regular expression
        expected_string_form_4 = r'^\( \+ \(1\)\*x0\) \*\* \( \+ \(\-1.5(|000000[0-9]*)\)\)$'
        self.assertRegexpMatches(str(poly_pow_float), expected_string_form_4)
        self.assertTrue(isinstance(poly_pow_float.exponent.coeffs[0], sp.Expr)) # ``-1`` should be converted to sympy expression

    #@attr('active')
    def test_pow_positive_integer(self):
        target_p0_zeroth_power = 0 * self.p0 + 1
        p0_zeroth_power = self.p0 ** 0
        self.assertEqual(sympify_expression(target_p0_zeroth_power - p0_zeroth_power), 0)

        target_p0_third_power = self.p0 * self.p0 * self.p0
        p0_third_power = self.p0 ** 3
        self.assertEqual(sympify_expression(target_p0_third_power - p0_third_power), 0)

        target_p1_sixth_power = self.p1 * self.p1 * self.p1 * self.p1 * self.p1 * self.p1
        p1_sixth_power = self.p1 ** 6
        self.assertEqual(sympify_expression(target_p1_sixth_power - p1_sixth_power), 0)

    def test_sympy_binding(self):
        a,b = sp.symbols('a b')

        p = Polynomial([(1,0),(0,1)],[a,b])
        p_squared = p * p
        self.assertEqual(str(p), ' + (a)*x0 + (b)*x1')
        self.assertEqual( (sympify_expression(str(p_squared)) - sympify_expression(' + (a**2)*x0**2 + (2*a*b)*x0*x1 + (b**2)*x1**2')).simplify() , 0)

        # should have minimal number of terms
        self.assertEqual(len(p_squared.coeffs), 3)
        self.assertEqual(p_squared.expolist.shape, (3,2))

        p_sum = p_squared + Polynomial([(1,1),(0,1)],[a,b])
        self.assertEqual( (sympify_expression(str(p_sum)) - sympify_expression(' + (a**2)*x0**2 + (2*a*b + a)*x0*x1 + (b**2)*x1**2 + (b)*x1')).simplify() , 0)

        # should have minimal number of terms
        self.assertEqual(len(p_sum.coeffs), 4)
        self.assertEqual(p_sum.expolist.shape, (4,2))


    def test_empty_expolist(self):
        polynomial = Polynomial([(0,1),(1,0),(2,1),(0,0)],[0,0,0,0])
        polynomial.simplify()
        self.assertGreater(len(polynomial.expolist), 0)
        self.assertGreater(len(polynomial.coeffs), 0)

        np.testing.assert_array_equal(polynomial.expolist, [[0,0]])
        np.testing.assert_array_equal(polynomial.coeffs, [0])

    def test_creation_from_expression(self):
        x,y, a,b,c = sp.symbols('x y a b c')
        polynomial_expression = a*x + b*y + c*x**2*y

        poly1 = Polynomial.from_expression(polynomial_expression, [x,y])
        poly2 = Polynomial.from_expression('a*x + b*y + c*x**2*y', ['x','y'])

        self.assertEqual((sympify_expression(str(poly1)) - polynomial_expression).simplify(), 0)
        self.assertEqual((sympify_expression(str(poly2)) - polynomial_expression).simplify(), 0)

        self.assertRaisesRegexp(TypeError, "\'x\*y\' is not.*symbol", Polynomial.from_expression, 'a*x + b*y + c*x**2*y', [x,x*y])
        self.assertRaisesRegexp(TypeError, "polysymbols.*at least one.*symbol", Polynomial.from_expression, 'a*x + b*y + c*x**2*y', [])

    def test_negative_exponent(self):
        vars = ['x','y']
        p1 = Polynomial.from_expression('A*y**-1*x + B*x*y*x**-2 + C*y**-1', vars)
        p2 = Polynomial([[1,-1],[-1,1],[0,-1]],['A','B','C'],vars)
        np.testing.assert_array_equal(p1.coeffs, p2.coeffs)
        np.testing.assert_array_equal(p1.expolist, p2.expolist)
        np.testing.assert_array_equal(p1.polysymbols, p2.polysymbols)

class TestExponentiatedPolynomial(unittest.TestCase):
    def test_init(self):
        ExponentiatedPolynomial([(1,2),(1,0),(2,1)], ['x',2,3])
        ExponentiatedPolynomial([(1,2),(1,0),(2,1)], ['x',2,3], exponent='A + eps')

    def test_string_form(self):
        polynomial1 = ExponentiatedPolynomial([(1,2),(1,0),(2,1)], ['x',2,3], polysymbols='y')
        string_polynomial1 = ' + (x)*y0*y1**2 + (2)*y0 + (3)*y0**2*y1'
        # if exponent is one, do not show it
        self.assertEqual(str(polynomial1), string_polynomial1)
        self.assertEqual(repr(polynomial1), string_polynomial1)

        polynomial2 = ExponentiatedPolynomial([(1,2),(1,0),(2,1)], ['x',2,3], polysymbols='y', exponent='A + eps')
        string_polynomial2 = '( + (x)*y0*y1**2 + (2)*y0 + (3)*y0**2*y1)**(A + eps)'
        self.assertEqual(str(polynomial2), string_polynomial2)
        self.assertEqual(repr(polynomial2), string_polynomial2)

    def test_copy(self):
        polynomial1 = ExponentiatedPolynomial([(0,1),(1,0),(2,1),(0,0)],['A','B','C','D'],exponent='eps')
        polynomial2 = polynomial1.copy()

        self.assertEqual(str(polynomial1), str(polynomial2))
        self.assertEqual(polynomial1.exponent, polynomial2.exponent)

    #@attr('active')
    def test_derive(self):
        A, B = sp.symbols('A B')
        polynomial = ExponentiatedPolynomial([(2,1),(0,0)],[A, B],exponent=sympify_expression('a + b*eps'))

        derivative_0 = sympify_expression( str(polynomial.derive(0)) )
        target_derivative_0 = sympify_expression('(a + b*eps)*(A*x0**2*x1 + B)**(a + b*eps - 1) * (2*A*x0*x1)')
        self.assertEqual( (derivative_0 - target_derivative_0).simplify() , 0 )

        derivative_1 = sympify_expression( str(polynomial.derive(1)) )
        target_derivative_1 = sympify_expression('(a + b*eps)*(A*x0**2*x1 + B)**(a + b*eps - 1) * (A*x0**2)')
        self.assertEqual( (derivative_1 - target_derivative_1).simplify() , 0 )


        polynomial = ExponentiatedPolynomial([(2,1),(0,0)],[A, B],exponent=Polynomial.from_expression('a + b*x0',['x0','x1']))
        derivative_0 = sympify_expression( str(polynomial.derive(0)) )
        target_derivative_0 = sympify_expression('(a + b*x0)*(A*x0**2*x1 + B)**(a + b*x0 - 1) * (2*A*x0*x1)   +   (A*x0**2*x1 + B)**(a + b*x0)*b*log(A*x0**2*x1 + B)')
        self.assertEqual( (derivative_0 - target_derivative_0).simplify() , 0 )

    #@attr('active')
    def test_derive_with_coeffs(self):
        expr = ExponentiatedPolynomial([(1,1),(0,0)],[Polynomial([(1,0)], ['A']), 'B'],exponent=Polynomial.from_expression('a + b*x0',['x0','x1']))
        derivative_0 = sympify_expression(expr.derive(0))
        target_derivative_0 = sympify_expression('(a + b*x0)*(A*x0**2*x1 + B)**(a + b*x0 - 1) * (2*A*x0*x1)   +   (A*x0**2*x1 + B)**(a + b*x0)*b*log(A*x0**2*x1 + B)')
        self.assertEqual( (derivative_0 - target_derivative_0).simplify() , 0 )

    #@attr('active')
    def test_simplify_zero_in_exponent(self):
        A, B = sp.symbols('A B')

        # <something>**0 = 1
        polynomial_to_power_zero_polynomial_exponent = ExponentiatedPolynomial([(2,1),(0,0)],[A, B],exponent=Polynomial.from_expression('0',['x0','x1'])).simplify()
        polynomial_to_power_zero_sympy_exponent = ExponentiatedPolynomial([(2,1),(0,0)],[A, B],exponent=sympify_expression('x-x')).simplify()
        polynomial_to_power_zero_numerical_exponent = ExponentiatedPolynomial([(2,1),(0,0)],[A, B],exponent=0).simplify()

        for p in (polynomial_to_power_zero_polynomial_exponent, polynomial_to_power_zero_sympy_exponent, polynomial_to_power_zero_numerical_exponent):
            self.assertTrue(type(p) is Polynomial)
            np.testing.assert_array_equal(p.coeffs, [1])
            np.testing.assert_array_equal(p.expolist, [[0,0]])

    #@attr('active')
    def test_simplify_one_in_exponent(self):
        A, B = sp.symbols('A B')

        # <something>**1 = <something>
        polynomial_to_power_one_polynomial_exponent = ExponentiatedPolynomial([(2,1),(0,0)],[A, B],exponent=Polynomial([[0,0,0],[0,0,0]], [2, -1])).simplify()
        polynomial_to_power_one_sympy_exponent = ExponentiatedPolynomial([(2,1),(0,0)],[A, B],exponent=sympify_expression('1/2+1/2')).simplify()
        polynomial_to_power_one_numerical_exponent = ExponentiatedPolynomial([(2,1),(0,0)],[A, B],exponent=1).simplify()

        for p in (polynomial_to_power_one_polynomial_exponent, polynomial_to_power_one_sympy_exponent, polynomial_to_power_one_numerical_exponent):
            self.assertTrue(type(p) is Polynomial)
            np.testing.assert_array_equal(p.coeffs, [A,B])
            np.testing.assert_array_equal(p.expolist, [(2,1),(0,0)])

    #@attr('active')
    def test_simplify_one_in_base(self):
        A = sp.symbols('A')

        # 1**<something> = 1
        polynomial_exponent = ExponentiatedPolynomial([(0,0),(0,0)],[A, 1-A],exponent=Polynomial([[0,9,0],[1,2,3]], [2, A])).simplify()
        sympy_exponent = ExponentiatedPolynomial([(0,0),(0,0)],[A, 1-A],exponent=sympify_expression('1/9+11/2 * eps')).simplify()
        numerical_exponent = ExponentiatedPolynomial([(0,0),(0,0)],[A, 1-A],exponent=np.pi).simplify()

        for p in (polynomial_exponent, sympy_exponent, numerical_exponent):
            self.assertTrue(type(p) is Polynomial)
            np.testing.assert_array_equal(p.coeffs, [1])
            np.testing.assert_array_equal(p.expolist, [(0,0)])

class TestProduct(unittest.TestCase):
    def test_init(self):
        p0 = Polynomial([(0,1),(1,0),(2,1)],['A','B','C'])
        p1 = Polynomial([(8,1),(1,5),(2,1)],['D','E','F'])
        p2 = Polynomial([(1,0,1),(2,1,0),(0,2,1)],['G','H','I'])

        # mismatch in number of parameters
        self.assertRaisesRegexp(TypeError, 'same number of variables.*all.*factors', Product,p0,p2)

        self.assertRaisesRegexp(AssertionError, 'at least one factor', Product)

        # Proper instantiation
        prod = Product(p0,p1)

        # made a copy?
        self.assertEqual(prod.factors[0].expolist[0,0],0)
        self.assertEqual(p0.expolist[0,0],0)
        prod.factors[0].expolist[0,0] = 5
        self.assertEqual(prod.factors[0].expolist[0,0],5)
        self.assertEqual(p0.expolist[0,0],0)

    def test_copy(self):
        p0 = Polynomial([(0,1),(1,0),(2,1)],['A','B','C'])
        p1 = Polynomial([(8,1),(1,5),(2,1)],['D','E','F'])

        orig = Product(p0,p1)
        copy = orig.copy()

        self.assertEqual(orig.factors[0].expolist[0,0],0)
        self.assertEqual(copy.factors[0].expolist[0,0],0)
        orig.factors[0].expolist[0,0] = 5
        self.assertEqual(orig.factors[0].expolist[0,0],5)
        self.assertEqual(copy.factors[0].expolist[0,0],0)

    #@attr('active')
    def test_derive(self):
        A, B = sp.symbols('A B')
        polynomial1 = Polynomial([(2,1),(0,0)],[A, B])
        polynomial2 = Polynomial([(0,0),(1,0)],[1, 2])
        prod = Product(polynomial1, polynomial2)

        derivative_0 = prod.derive(0)
        derivative_0_1 = sympify_expression( str(derivative_0.derive(1)) )
        target_derivative_0_1 = sympify_expression('(A*x0**2*x1 + B) * (1 + 2*x0)')
        target_derivative_0_1 = sympify_expression('(2*A*x0*x1) * (1 + 2*x0) + 2*(A*x0**2*x1 + B)')
        target_derivative_0_1 = sympify_expression('(2*A*x0) * (1 + 2*x0) + 2*(A*x0**2)')
        self.assertEqual( (derivative_0_1 - target_derivative_0_1).simplify() , 0 )

    def test_string_form(self):
        p0 = ExponentiatedPolynomial([(0,1)],['A'],exponent='exponent')
        p1 = Polynomial([(8,1),(1,5),(2,1)],['B','C','D'])
        prod = Product(p0,p1)
        string_prod = '(( + (A)*x1)**(exponent)) * ( + (B)*x0**8*x1 + (C)*x0*x1**5 + (D)*x0**2*x1)'

        self.assertEqual(str(prod), string_prod)
        self.assertEqual(repr(prod), string_prod)

    #@attr('active')
    def test_simplify(self):
        complicated_one = Polynomial([(2,1),(0,0),(2,1),(0,0)],['A','B','-A','-B+1'])
        prod = Product(complicated_one,complicated_one,complicated_one)

        simplified_prod = prod.simplify()

        self.assertEqual( sympify_expression(simplified_prod) , 1 )
        self.assertEqual( sympify_expression(prod) , 1 )
        self.assertTrue(type(prod) is Product)
        self.assertEqual(len(prod.factors), 1)

class TestProductRule(unittest.TestCase):
    def setUp(self):
        self.poly1 = Polynomial.from_expression('x+x*y', ['x','y'])
        self.poly2 = Polynomial.from_expression('x**2-y**2', ['x','y'])
        self.p0 = ProductRule(self.poly1, self.poly2)
        self.target_dp0_dx = sympify_expression('(1+y) * (x**2-y**2) + (x+x*y) * (2*x)')
        self.target_ddp0_dxdy = sympify_expression('(x**2-y**2) + (1+y) * (-2*y)   +   x * (2*x)')
        self.target_dddp0_dxdydx = sympify_expression('2*x + 4*x')

    #@attr('active')
    def test_string_form_basic(self):
        target_str_p0 = ' + (1) * ( + (1)*x + (1)*x*y) * ( + (-1)*y**2 + (1)*x**2)'
        self.assertEqual(str(self.p0), target_str_p0)
        self.assertEqual(repr(self.p0), target_str_p0)

    #@attr('active')
    def test_copy(self):
        p0 = ProductRule(self.poly1, self.poly2)
        p1 = p0.copy()

        self.assertEqual(len(p0.expressions), len(p1.expressions))

        for expr0, expr1 in zip(p0.expressions, p1.expressions):
            for key in list(expr0.keys()) + list(expr1.keys()):
                self.assertFalse(expr0[key] is expr1[key])

    #@attr('active')
    def test_derive(self):
        dp0_dx = self.p0.derive(0)
        self.assertEqual(  (sympify_expression(dp0_dx) - self.target_dp0_dx).simplify()  ,   0   )

        ddp0_dxdy = self.p0.derive(0).derive(1)
        self.assertEqual(  (sympify_expression(ddp0_dxdy) - self.target_ddp0_dxdy).simplify()  ,   0   )

        dddp0_dxdydx = self.p0.derive(0).derive(1).derive(0)
        self.assertEqual(  (sympify_expression(dddp0_dxdydx) - self.target_dddp0_dxdydx).simplify()  ,   0   )

    #@attr('active')
    def test_simplify(self):
        ddp0_dx_dx = self.p0.derive(0).derive(0)
        simplified_ddp0_dx_dx = self.p0.derive(0).derive(0).simplify()
        target_ddp0_dx_dx = sympify_expression('(1+y) * (2*x) + 4 * x * (1+y)')

        self.assertEqual(   (sympify_expression(simplified_ddp0_dx_dx) - target_ddp0_dx_dx).simplify()   ,   0   )
        self.assertLess(len(simplified_ddp0_dx_dx.coeffs), len(ddp0_dx_dx.coeffs))

    #@attr('active')
    def test_simplify_to_zero(self):
        p0 = Polynomial([[1,2,0]], [1])
        p1 = Polynomial([[2,1,0]], [1])

        prod_rule = ProductRule(p0,p1).derive(-1)
        prod_rule_simplify_returned = prod_rule.simplify()

        self.assertFalse(prod_rule is prod_rule_simplify_returned)
        self.assertTrue(type(prod_rule) is ProductRule)
        self.assertTrue(type(prod_rule_simplify_returned) is Polynomial)
        self.assertEqual( sympify_expression(prod_rule).simplify(), 0 )
        self.assertEqual( sympify_expression(prod_rule_simplify_returned).simplify(), 0 )

    #@attr('active')
    def test_to_sum(self):
        p0_sum = self.p0.copy().to_sum()
        self.assertTrue(type(p0_sum) is Sum)
        self.assertEqual(   (sympify_expression(p0_sum) - sympify_expression(self.p0)).simplify()   ,   0   )

    #@attr('active')
    def test_replace(self):
        z = sp.symbols('z')
        dddp0_dxdydx = self.p0.derive(0).derive(1).derive(0).simplify()
        replaced = dddp0_dxdydx.replace(0,z)
        removed = dddp0_dxdydx.replace(0,z, True)

        self.assertEqual(replaced.symbols, sp.symbols(['x','y']))
        self.assertEqual(removed.symbols, sp.symbols(['y']))

        for derivative in [removed, replaced]:
            self.assertEqual(   (sympify_expression(derivative) - self.target_dddp0_dxdydx.subs('x','z')).simplify()   ,   0   )

#@attr('active')
class TestPow(unittest.TestCase):
    def setUp(self):
        self.zero = Polynomial.from_expression(0, ['x0','x1','x2'])
        self.one = Polynomial.from_expression(1, ['x0','x1','x2'])
        self.something = Polynomial([(0,0,1),(0,1,0),(1,0,0)], ['a','b','c'])

    def test_init(self):
        p0 = Polynomial([(0,1),(1,0),(2,1)],['A','B','C'])
        p1 = Polynomial([(8,1),(1,5),(2,1)],['D','E','F'])
        p2 = Polynomial([(1,0,1),(2,1,0),(0,2,1)],['G','H','I'])

        # mismatch in number of parameters
        self.assertRaisesRegexp(TypeError, 'same number of variables.*base.*exponent', Pow,p0,p2)

        # Proper instantiation
        exp = Pow(p0,p1)

        # made a copy?
        self.assertEqual(exp.base.expolist[0,0],0)
        self.assertEqual(p0.expolist[0,0],0)
        exp.base.expolist[0,0] = 5
        self.assertEqual(exp.base.expolist[0,0],5)
        self.assertEqual(p0.expolist[0,0],0)

        self.assertEqual(exp.exponent.expolist[0,0],8)
        self.assertEqual(p1.expolist[0,0],8)
        exp.exponent.expolist[0,0] = 5
        self.assertEqual(exp.exponent.expolist[0,0],5)
        self.assertEqual(p1.expolist[0,0],8)

    def test_copy(self):
        p0 = Polynomial([(0,1),(1,0),(2,1)],['A','B','C'])
        p1 = Polynomial([(8,1),(1,5),(2,1)],['D','E','F'])

        orig = Pow(p0,p1)
        copy = orig.copy()

        self.assertEqual(orig.base.expolist[0,0],0)
        self.assertEqual(copy.base.expolist[0,0],0)
        orig.base.expolist[0,0] = 5
        self.assertEqual(orig.base.expolist[0,0],5)
        self.assertEqual(copy.base.expolist[0,0],0)

        self.assertEqual(orig.exponent.expolist[0,0],8)
        self.assertEqual(copy.exponent.expolist[0,0],8)
        orig.exponent.expolist[0,0] = 5
        self.assertEqual(orig.exponent.expolist[0,0],5)
        self.assertEqual(copy.exponent.expolist[0,0],8)

    #@attr('active')
    def test_simplify_zero_in_exponent(self):
        # <something>**0 = 1
        one = Pow(self.something, self.zero).simplify()
        self.assertTrue(type(one) is Polynomial)
        np.testing.assert_array_equal(one.coeffs, self.one.coeffs)
        np.testing.assert_array_equal(one.expolist, self.one.expolist)

    #@attr('active')
    def test_simplify_one_in_exponent(self):
        # <something>**1 = <something>
        poly = Pow(self.something, self.one).simplify()
        self.assertTrue(type(poly) is Polynomial)
        np.testing.assert_array_equal(poly.coeffs, self.something.coeffs)
        np.testing.assert_array_equal(poly.expolist, self.something.expolist)

    #@attr('active')
    def test_simplify_one_in_base(self):
        # 1**<something> = 1
        complicated_one = Polynomial([[1,2,3],[1,2,3],[0,0,0],[0,0,0]], ['A','-A','B','1-B'])
        complicated_one = Product(complicated_one, complicated_one, self.one) ** sp.symbols('exponent1')
        one = Pow(complicated_one, self.something**sp.symbols('exponent2')).simplify()
        self.assertTrue(type(one) is Polynomial)
        np.testing.assert_array_equal(one.coeffs, self.one.coeffs)
        np.testing.assert_array_equal(one.expolist, self.one.expolist)

    #@attr('active')
    def test_simplify_remove_expression(self):
        constant_poly = Polynomial([[0,0,0], [0,0,0], [0,0,0]], ['a*3', 4*9, 'b'])
        base = Polynomial([(0,0,1),(0,1,0),(1,0,0)], ['a','b','c'])

        exponentiated = Pow(base, constant_poly).simplify()

        self.assertTrue(type(exponentiated) is ExponentiatedPolynomial)
        self.assertFalse( isinstance(exponentiated.exponent, _Expression) )
        self.assertEqual( (sympify_expression(exponentiated) - sympify_expression('(a*x2 + b*x1 + c*x0)**(a*3 + 4*9 + b)')).simplify() , 0 )

    #@attr('active')
    def test_derive(self):
        A, B = sp.symbols('A B')
        polynomial1 = Polynomial.from_expression('A*x0 + B*x1', ['x0','x1'])
        polynomial2 = Polynomial.from_expression('x1', ['x0','x1'])
        exp = Pow(polynomial1, polynomial2)

        derivative_0 = exp.derive(0)
        target_derivative_0 = sympify_expression('(A*x0 + B*x1) ** (x1 - 1) * A*x1')
        self.assertEqual( (sympify_expression(derivative_0) - target_derivative_0).simplify() , 0 )

        derivative_0_1 = sympify_expression(derivative_0.derive(1))
        target_derivative_0_1 = sympify_expression('A*x1 *( (A*x0 + B*x1) ** (x1 - 1) * (log(A*x0 + B*x1) + (x1 - 1)*B/(A*x0 + B*x1)) ) + (A*x0 + B*x1) ** (x1 - 1) * A')
        self.assertEqual( (derivative_0_1 - target_derivative_0_1).simplify() , 0 )

        # `_Expression` in exponent
        exp = Pow(polynomial1, Sum(polynomial1, polynomial2))
        derivative_0 = exp.derive(0)
        target_derivative_0 = sympify_expression('(A*x0 + B*x1)**(A*x0 + B*x1 + x1) * ( A*log(A*x0 + B*x1) + (A*x0 + B*x1 + x1)/(A*x0 + B*x1)*A )')
        self.assertEqual( (sympify_expression(derivative_0) - target_derivative_0).simplify() , 0 )

    def test_string_form(self):
        p0 = ExponentiatedPolynomial([(0,1)],['A'],exponent='exponent')
        p1 = Polynomial([(8,1),(1,5),(2,1)],['B','C','D'])
        exp = Pow(p0,p1)
        string_pow = '(( + (A)*x1)**(exponent)) ** ( + (B)*x0**8*x1 + (C)*x0*x1**5 + (D)*x0**2*x1)'

        self.assertEqual(str(exp), string_pow)
        self.assertEqual(repr(exp), string_pow)

#@attr('active')
class TestLog(unittest.TestCase):
    def test_init(self):
        p0 = Polynomial([(0,1),(1,0),(2,1)],['A','B','C'])

        # Proper instantiation
        ln = Log(p0)

        # made a copy?
        self.assertEqual(ln.arg.expolist[0,0],0)
        self.assertEqual(p0.expolist[0,0],0)
        ln.arg.expolist[0,0] = 5
        self.assertEqual(ln.arg.expolist[0,0],5)
        self.assertEqual(p0.expolist[0,0],0)

    def test_copy(self):
        p0 = Polynomial([(0,1),(1,0),(2,1)],['A','B','C'])

        orig = Log(p0)
        copy = orig.copy()

        self.assertEqual(orig.arg.expolist[0,0],0)
        self.assertEqual(copy.arg.expolist[0,0],0)
        orig.arg.expolist[0,0] = 5
        self.assertEqual(orig.arg.expolist[0,0],5)
        self.assertEqual(copy.arg.expolist[0,0],0)

    #@attr('active')
    def test_simplify(self):
        one = Polynomial.from_expression(1, ['x0','x1','x2'])

        zero = Log(one).simplify()

        self.assertTrue(type(zero) is Polynomial)
        self.assertEqual(sympify_expression(zero), 0)
        np.testing.assert_array_equal(zero.coeffs, [0])
        np.testing.assert_array_equal(zero.expolist, [[0,0,0]])

    #@attr('active')
    def test_simplify_log_of_polynomial_one(self):
        one = Polynomial([[4,1,5],[0,0,0]], [0,Polynomial.from_expression(1, ['x','y','z'])], ['x','y','z'])

        zero = Log(one).simplify()

        self.assertTrue(type(zero) is Polynomial)
        self.assertEqual(sympify_expression(zero), 0)
        np.testing.assert_array_equal(zero.coeffs, [0])
        np.testing.assert_array_equal(zero.expolist, [[0,0,0]])

    def test_derive(self):
        polynomial = Polynomial.from_expression('A*x0 + B*x1', ['x0','x1'])
        ln = Log(polynomial)

        derivative_0 = ln.derive(0).simplify()
        target_derivative_0 = sympify_expression('1/(A*x0 + B*x1)*A')
        self.assertEqual( (sympify_expression(derivative_0) - target_derivative_0).simplify() , 0 )

        derivative_0_1 = sympify_expression(derivative_0.derive(1))
        target_derivative_0_1 = sympify_expression('-A * (A*x0 + B*x1)**(-2) * B')
        self.assertEqual( (derivative_0_1 - target_derivative_0_1).simplify() , 0 )

    def test_string_form(self):
        p1 = Polynomial([(8,1),(1,5),(2,1)],['B','C','D'])
        ln = Log(p1)
        string_ln = 'log( + (B)*x0**8*x1 + (C)*x0*x1**5 + (D)*x0**2*x1)'

        self.assertEqual(str(ln), string_ln)
        self.assertEqual(repr(ln), string_ln)

class TestSum(unittest.TestCase):
    def test_init(self):
        p0 = Polynomial([(0,1),(1,0),(2,1)],['A','B','C'])
        p1 = Polynomial([(8,1),(1,5),(2,1)],['D','E','F'])
        p2 = Polynomial([(1,0,1),(2,1,0),(0,2,1)],['G','H','I'])

        # mismatch in number of parameters
        self.assertRaisesRegexp(TypeError, 'same number of variables.*all.*summands', Sum,p0,p2)

        self.assertRaisesRegexp(AssertionError, 'at least one summand', Sum)

        # Proper instantiation
        psum = Sum(p0,p1)

        # made a copy?
        self.assertEqual(psum.summands[0].expolist[0,0],0)
        self.assertEqual(p0.expolist[0,0],0)
        psum.summands[0].expolist[0,0] = 5
        self.assertEqual(psum.summands[0].expolist[0,0],5)
        self.assertEqual(p0.expolist[0,0],0)

    def test_copy(self):
        p0 = Polynomial([(0,1),(1,0),(2,1)],['A','B','C'])
        p1 = Polynomial([(8,1),(1,5),(2,1)],['D','E','F'])

        orig = Sum(p0,p1)
        copy = orig.copy()

        self.assertEqual(orig.summands[0].expolist[0,0],0)
        self.assertEqual(copy.summands[0].expolist[0,0],0)
        orig.summands[0].expolist[0,0] = 5
        self.assertEqual(orig.summands[0].expolist[0,0],5)
        self.assertEqual(copy.summands[0].expolist[0,0],0)

    def test_string_form(self):
        p0 = ExponentiatedPolynomial([(0,1)],['A'],exponent='exponent')
        p1 = Polynomial([(8,1),(1,5),(2,1)],['B','C','D'])
        sum = Sum(p0,p1)
        string_sum = '(( + (A)*x1)**(exponent)) + ( + (B)*x0**8*x1 + (C)*x0*x1**5 + (D)*x0**2*x1)'

        self.assertEqual(str(sum), string_sum)
        self.assertEqual(repr(sum), string_sum)

    def test_derive(self):
        A, B = sp.symbols('A B')

        p0 = ExponentiatedPolynomial([(0,1)],[A])
        p1 = ExponentiatedPolynomial([(2,1)],[B])
        psum = Sum(p0,p1)

        derivative_0 = sympify_expression( psum.derive(0) )
        target_derivative_0 = sympify_expression( '2*B*x0*x1' )
        self.assertEqual( (derivative_0 - target_derivative_0).simplify() , 0 )

    #@attr('active')
    def test_simplify(self):
        complicated_zero = Polynomial([(2,1),(0,0),(2,1),(0,0)],['A','B-1','-A','-B+1'])
        psum = Sum(complicated_zero,complicated_zero,complicated_zero)

        simplified_psum = psum.simplify()

        self.assertEqual( sympify_expression(simplified_psum) , 0 )
        self.assertEqual( sympify_expression(psum) , 0 )
        self.assertTrue(type(psum) is Sum)
        self.assertEqual(len(psum.summands), 1)

class TestLogOfPolynomial(unittest.TestCase):
    def test_string_form(self):
        p0 = LogOfPolynomial([(0,1),(1,0),(2,1)],['A','B','C'])
        str_p0 = 'log( + (A)*x1 + (B)*x0 + (C)*x0**2*x1)'

        self.assertEqual(str(p0),str_p0)
        self.assertEqual(repr(p0),str_p0)

    def test_construct_from_expression(self):
        p1 = LogOfPolynomial.from_expression('D*x0**8*x1 + E*x0*x1**5 + F*x0*x0*x1',['x0','x1'])
        str_p1 = 'log( + (D)*x0**8*x1 + (E)*x0*x1**5 + (F)*x0**2*x1)'
        sympy_p1 = sympify_expression(str_p1)

        self.assertEqual(str(p1),repr(p1))
        self.assertEqual( sympify_expression(repr(p1)) - sympy_p1 , 0 )

    def test_derive(self):
        expr = LogOfPolynomial([(2,1),(0,1)],['A', 'B'])

        derivative_0 = expr.derive(0)
        sympified_derivative_0 = sympify_expression( str(derivative_0) )
        target_derivative_0 = sympify_expression('1/(A*x0**2*x1 + B*x1) * 2*A*x0*x1')
        self.assertEqual( (sympified_derivative_0 - target_derivative_0).simplify() , 0 )
        self.assertEqual(type(derivative_0), Product)
        self.assertEqual(len(derivative_0.factors), 2)
        self.assertEqual(type(derivative_0.factors[0]), ExponentiatedPolynomial)

        derivative_1 = sympify_expression( str(expr.derive(1)) )
        target_derivative_1 = sympify_expression('1/(A*x0**2*x1 + B*x1) * (A*x0**2 + B)')
        self.assertEqual( (derivative_1 - target_derivative_1).simplify() , 0 )

    #@attr('active')
    def test_derive_with_coeffs(self):
        expr = LogOfPolynomial([(1,1),(0,1)],[Polynomial([(1,0)], ['A']), 'B'])
        derivative_0 = sympify_expression(expr.derive(0))
        target_derivative_0 = sympify_expression('1/(A*x0**2*x1 + B*x1) * 2*A*x0*x1')
        self.assertEqual( (derivative_0 - target_derivative_0).simplify() , 0)

    def test_simplify(self):
        # log(1) = 0
        expr = LogOfPolynomial([(0,0),(0,0)],[1, 0], polysymbols=['x','y'])
        simplified_expr = expr.simplify()

        self.assertTrue(type(simplified_expr) is Polynomial)
        np.testing.assert_array_equal(simplified_expr.expolist, [[0,0]])
        np.testing.assert_array_equal(simplified_expr.coeffs, [0])

    #@attr('active')
    def test_simplify_log_of_polynomial_one(self):
        one = Polynomial.from_expression(1, ['x','y','z'])

        zero = LogOfPolynomial([[4,1,5],[0,0,0]], [0,one], ['x','y','z']).simplify()

        self.assertTrue(type(zero) is Polynomial)
        self.assertEqual(sympify_expression(zero), 0)
        np.testing.assert_array_equal(zero.coeffs, [0])
        np.testing.assert_array_equal(zero.expolist, [[0,0,0]])

#@attr('active')
class TestInsertion(unittest.TestCase):
    def test_insert_value_polynomial(self):
        poly = Polynomial([(0,0),(1,0),(0,1),(0,2)],['A','B','C','D'])
        replaced_poly = poly.replace(index=1,value=sympify_expression('1/2'))
        self.assertEqual( sympify_expression(str(replaced_poly)) - sympify_expression('A + B*x0 + C/2 + D/4') , 0 )
        self.assertEqual(replaced_poly.number_of_variables, 2)

        removed_poly = poly.replace(index=1,value=sympify_expression('1/2'), remove=True)
        self.assertEqual( sympify_expression(str(replaced_poly)) - sympify_expression('A + B*x0 + C/2 + D/4') , 0 )
        self.assertEqual(removed_poly.number_of_variables, 1)
        self.assertEqual(removed_poly.expolist.shape[1], 1)

    #@attr('active')
    def test_insert_polynomial_coeffs(self):
        poly = Polynomial([(0,0),(1,0),(0,1),(0,2)],['A','B',Polynomial([(0,1)], ['C']),'D'])

        replaced_poly = poly.replace(index=1,value=sympify_expression('1/2'))
        self.assertEqual( sympify_expression(str(replaced_poly)) - sympify_expression('A + B*x0 + C/4 + D/4') , 0 )
        self.assertEqual(replaced_poly.number_of_variables, 2)

        removed_poly = poly.replace(index=1,value=sympify_expression('1/2'), remove=True)
        self.assertEqual( sympify_expression(str(replaced_poly)) - sympify_expression('A + B*x0 + C/4 + D/4') , 0 )
        self.assertEqual(removed_poly.number_of_variables, 1)
        self.assertEqual(removed_poly.expolist.shape[1], 1)

    #@attr('active')
    def test_insert_value_polynomial_negative_index(self):
        poly = Polynomial([(0,0),(1,0),(0,1),(0,2)],['A','B','C','D'])
        replaced_poly = poly.replace(index=-1,value=sympify_expression('1/2'))
        self.assertEqual( sympify_expression(str(replaced_poly)) - sympify_expression('A + B*x0 + C/2 + D/4') , 0 )
        self.assertEqual(replaced_poly.number_of_variables, 2)
        self.assertEqual(replaced_poly.polysymbols, sympify_expression(['x0','x1']))

        removed_poly = poly.replace(index=-1,value=sympify_expression('1/2'), remove=True)
        self.assertEqual( sympify_expression(str(replaced_poly)) - sympify_expression('A + B*x0 + C/2 + D/4') , 0 )
        self.assertEqual(removed_poly.number_of_variables, 1)
        self.assertEqual(removed_poly.expolist.shape[1], 1)
        self.assertEqual(removed_poly.polysymbols, sympify_expression(['x0']))

    def test_insert_value_polynomial_product(self):
        poly0 = Polynomial([(0,0),(1,0),(0,1),(0,2)],['A','B','C','D'])
        poly1 = Polynomial([(0,0),(5,0)],['E',1])
        prod = Product(poly0,poly1)
        replaced_prod = prod.replace(index=1,value=0)
        self.assertEqual( (sympify_expression(str(replaced_prod)) - sympify_expression('(A + B*x0) * (E + x0**5)')).simplify() , 0 )
        self.assertEqual(replaced_prod.number_of_variables, 2)

        removed_prod = prod.replace(index=1, value=0, remove=True)
        self.assertEqual( (sympify_expression(str(removed_prod)) - sympify_expression('(A + B*x0) * (E + x0**5)')).simplify() , 0 )
        self.assertEqual(removed_prod.number_of_variables, 1)
        self.assertEqual(removed_prod.factors[0].expolist.shape[1], 1)
        self.assertEqual(removed_prod.factors[1].expolist.shape[1], 1)

    def test_insert_value_polynomial_sum(self):
        poly0 = Polynomial([(0,0),(1,0),(0,1),(0,2)],['A','B','C','D'])
        poly1 = Polynomial([(0,0),(5,0)],['E',1])
        polysum = Sum(poly0,poly1)
        replaced_sum = polysum.replace(index=1,value=0)
        self.assertEqual( (sympify_expression(str(replaced_sum)) - sympify_expression('A + B*x0 + E + x0**5')).simplify() , 0 )
        self.assertEqual(replaced_sum.number_of_variables, 2)

        removed_sum = polysum.replace(index=1,value=0, remove=True)
        self.assertEqual( (sympify_expression(str(removed_sum)) - sympify_expression('(A + B*x0) + (E + x0**5)')).simplify() , 0 )
        self.assertEqual(removed_sum.number_of_variables, 1)
        self.assertEqual(removed_sum.summands[0].expolist.shape[1], 1)
        self.assertEqual(removed_sum.summands[1].expolist.shape[1], 1)

    def test_insert_in_exponent(self):
        exponent = Polynomial([(0,0),(5,0)],['E',1])
        poly1 = ExponentiatedPolynomial([(0,0),(1,0),(0,1),(0,2)],['A','B','C','D'],exponent)
        replaced = poly1.replace(index=0,value=0)
        self.assertEqual( (sympify_expression(str(replaced)) - sympify_expression('(A + C*x1 + D*x1**2)**(E)')).simplify() , 0 )
        self.assertEqual(replaced.number_of_variables, 2)
        self.assertEqual(replaced.expolist.shape[1], 2)
        self.assertEqual(replaced.exponent.number_of_variables, 2)
        self.assertEqual(replaced.exponent.expolist.shape[1], 2)

        removed = poly1.replace(index=0, value=2, remove=True)
        self.assertEqual( (sympify_expression(str(removed)) - sympify_expression('(A + 2*B + C*x1 + D*x1**2)**(E + 2**5)')).simplify() , 0 )
        self.assertEqual(removed.number_of_variables, 1)
        self.assertEqual(removed.exponent.number_of_variables, 1)
        self.assertEqual(removed.expolist.shape[1], 1)
        self.assertEqual(removed.exponent.expolist.shape[1], 1)

    def test_insert_in_pow(self):
        exponent = Polynomial([(0,0),(5,0)],['E',1])
        base = Polynomial([(0,0),(1,0),(0,1),(0,2)],['A','B','C','D'])
        expr = Pow(base, exponent)

        replaced = expr.replace(index=0,value=0)
        self.assertEqual( (sympify_expression(str(replaced)) - sympify_expression('(A + C*x1 + D*x1**2)**(E)')).simplify() , 0 )
        self.assertEqual(replaced.number_of_variables, 2)

        removed = expr.replace(index=0, value=2, remove=True)
        self.assertEqual( (sympify_expression(str(removed)) - sympify_expression('(A + 2*B + C*x1 + D*x1**2)**(E + 2**5)')).simplify() , 0 )
        self.assertEqual(removed.number_of_variables, 1)

    #@attr('active')
    def test_insert_in_log(self):
        arg = Polynomial([(0,0),(1,0),(0,1),(0,2)],['A','B','C','D'])
        expr = Log(arg)

        replaced = expr.replace(index=1,value=sympify_expression('(3*k+1)'))
        self.assertEqual( (sympify_expression(str(replaced)) - sympify_expression('log(A + B*x0 + C*(3*k+1) + D*(3*k+1)**2)')).simplify() , 0 )
        self.assertEqual(replaced.number_of_variables, 2)

        removed = expr.replace(index=1, value=sympify_expression('(3*k+1)'), remove=True)
        self.assertEqual( (sympify_expression(str(removed)) - sympify_expression('log(A + B*x0 + C*(3*k+1) + D*(3*k+1)**2)')).simplify() , 0 )
        self.assertEqual(removed.number_of_variables, 1)

    #@attr('active')
    def test_insert_in_Function(self):
        exponent = Polynomial([(0,0),(5,0)],['E',1])
        poly1 = ExponentiatedPolynomial([(0,0),(1,0),(0,1),(0,2)],['A','B','C','D'],exponent)
        poly2 = Polynomial([(1,0),(0,1)],[1,1])
        f = Function('f', poly1, poly2)

        replaced = f.replace(index=0,value=0)
        self.assertEqual( (sympify_expression(str(replaced)) - sympify_expression('f( (A + C*x1 + D*x1**2)**(E), x1)')).simplify() , 0 )
        self.assertEqual(replaced.number_of_variables, 2)
        self.assertEqual(replaced.arguments[0].expolist.shape[1], 2)
        self.assertEqual(replaced.arguments[1].expolist.shape[1], 2)
        self.assertEqual(replaced.arguments[1].number_of_variables, 2)
        self.assertEqual(replaced.arguments[0].exponent.number_of_variables, 2)
        self.assertEqual(replaced.arguments[0].exponent.expolist.shape[1], 2)

        removed = f.replace(index=0, value=2, remove=True)
        self.assertEqual( (sympify_expression(str(removed)) - sympify_expression('f( (A + 2*B + C*x1 + D*x1**2)**(E + 2**5), 2 + x1)')).simplify() , 0 )
        self.assertEqual(removed.number_of_variables, 1)
        self.assertEqual(removed.arguments[0].expolist.shape[1], 1)
        self.assertEqual(removed.arguments[1].expolist.shape[1], 1)
        self.assertEqual(removed.arguments[1].number_of_variables, 1)
        self.assertEqual(removed.arguments[0].exponent.number_of_variables, 1)
        self.assertEqual(removed.arguments[0].exponent.expolist.shape[1], 1)

#@attr('active')
class TestGetSymbols(unittest.TestCase):
    def setUp(self):
        self.polysymbols = sympify_expression(['x','y'])
        self.p0 = Polynomial.from_expression('a*x + b*y', self.polysymbols)
        self.p1 = Polynomial.from_expression('c*x + d*y', self.polysymbols)

    def test_get_from_polynomial(self):
        self.assertEqual(self.p0.symbols, self.polysymbols)
        self.assertEqual(self.p1.symbols, self.polysymbols)

    def test_get_from_sum(self):
        expr = Sum(self.p0, self.p1)
        self.assertEqual(expr.symbols, self.polysymbols)

    def test_get_from_product(self):
        expr = Product(self.p0, self.p1)
        self.assertEqual(expr.symbols, self.polysymbols)

    def test_get_from_log(self):
        expr = Sum(Log(self.p0), Log(self.p1))
        self.assertEqual(expr.symbols, self.polysymbols)

    def test_get_from_pow(self):
        expr = Pow(Log(self.p0), Log(self.p1))
        self.assertEqual(expr.symbols, self.polysymbols)

    def test_get_from_function(self):
        expr = Function('f', self.p0, self.p1)
        self.assertEqual(expr.symbols, self.polysymbols)

class TestExpressionOperators(unittest.TestCase):
    def setUp(self):
        self.polysymbols = sympify_expression(['x','y'])
        self.p0 = ExponentiatedPolynomial([(1,0),(0,1)], ['a','b'], 'exponent', self.polysymbols)
        self.p1 = LogOfPolynomial([(1,0),(0,1)], ['c','d'], self.polysymbols)

    #@attr('active')
    def test_add(self):
        sum_p0_p1 = self.p0 + self.p1
        sympified_sum_p0_p1 = sympify_expression(sum_p0_p1)
        target_sum_p0_p1 = sympify_expression('(a*x + b*y)**exponent + log(c*x + d*y)')
        self.assertTrue(type(sum_p0_p1) is Sum)
        self.assertEqual( (sympified_sum_p0_p1 - target_sum_p0_p1).simplify() , 0 )

        sum_p1_p0 = self.p1 + self.p0
        sympified_sum_p1_p0 = sympify_expression(sum_p1_p0)
        target_sum_p1_p0 = target_sum_p0_p1
        self.assertTrue(type(sum_p1_p0) is Sum)
        self.assertEqual( (sympified_sum_p1_p0 - target_sum_p1_p0).simplify() , 0 )

        for expr, sympified_expr in zip([self.p0, self.p1], sympify_expression(['(a*x + b*y)**exponent', 'log(c*x + d*y)'])):
            # other operand is number
            sum_p0_one = sympify_expression(expr + 1)
            self.assertEqual( (sympified_expr + 1 - sum_p0_one).simplify() , 0)

            sum_one_p0 = sympify_expression(1 + expr)
            self.assertEqual( (sympified_expr + 1 - sum_one_p0).simplify() , 0)

            # other operand is sympy symbol
            sum_p0_K = sympify_expression(expr + sympify_expression('K'))
            self.assertEqual( (sympified_expr + sympify_expression('K') - sum_p0_K).simplify() , 0)

            sum_K_p0 = sympify_expression(sympify_expression('K') + expr)
            self.assertEqual( (sympified_expr + sympify_expression('K') - sum_K_p0).simplify() , 0)

    #@attr('active')
    def test_neg(self):
        minus_p0 = sympify_expression(-self.p0)
        target_minus_p0 = sympify_expression('-(a*x + b*y)**exponent')
        self.assertEqual( (minus_p0 - target_minus_p0).simplify() , 0)

        minus_p1 = sympify_expression(-self.p1)
        target_minus_p1 = sympify_expression('-log(c*x + d*y)')
        self.assertEqual( (minus_p1 - target_minus_p1).simplify() , 0)

    #@attr('active')
    def test_sub(self):
        diff_p0_p1 = self.p0 - self.p1
        sympified_diff_p0_p1 = sympify_expression(diff_p0_p1)
        target_diff_p0_p1 = sympify_expression('(a*x + b*y)**exponent - log(c*x + d*y)')
        self.assertTrue(type(diff_p0_p1) is Sum)
        self.assertEqual( (sympified_diff_p0_p1 - target_diff_p0_p1).simplify() , 0 )

        diff_p1_p0 = self.p1 - self.p0
        sympified_diff_p1_p0 = sympify_expression(diff_p1_p0)
        target_diff_p1_p0 = -target_diff_p0_p1
        self.assertTrue(type(diff_p1_p0) is Sum)
        self.assertEqual( (sympified_diff_p1_p0 - target_diff_p1_p0).simplify() , 0 )

        for expr, sympified_expr in zip([self.p0, self.p1], sympify_expression(['(a*x + b*y)**exponent', 'log(c*x + d*y)'])):
            # other operand is number
            diff_p0_one = sympify_expression(expr -1)
            self.assertEqual( (sympified_expr - 1 - diff_p0_one).simplify() , 0)

            diff_one_p0 = sympify_expression(1 - expr)
            self.assertEqual( (1 - sympified_expr - diff_one_p0).simplify() , 0)

            # other operand is sympy symbol
            diff_p0_K = sympify_expression(expr - sympify_expression('K'))
            self.assertEqual( (sympified_expr - sympify_expression('K') - diff_p0_K).simplify() , 0)

            diff_K_p0 = sympify_expression(sympify_expression('K') - expr)
            self.assertEqual( (sympify_expression('K') - sympified_expr - diff_K_p0).simplify() , 0)

    #@attr('active')
    def test_mul(self):
        prod_p0_p1 = self.p0 * self.p1
        sympified_prod_p0_p1 = sympify_expression(prod_p0_p1)
        target_prod_p0_p1 = sympify_expression('(a*x + b*y)**exponent * log(c*x + d*y)')
        self.assertTrue(type(prod_p0_p1) is Product)
        self.assertEqual( (sympified_prod_p0_p1 - target_prod_p0_p1).simplify() , 0 )

        prod_p1_p0 = self.p1 * self.p0
        sympified_prod_p1_p0 = sympify_expression(prod_p1_p0)
        target_prod_p1_p0 = target_prod_p0_p1
        self.assertTrue(type(prod_p1_p0) is Product)
        self.assertEqual( (sympified_prod_p1_p0 - target_prod_p1_p0).simplify() , 0 )

        for expr, sympified_expr in zip([self.p0, self.p1], sympify_expression(['(a*x + b*y)**exponent', 'log(c*x + d*y)'])):
            # other operand is number
            prod_p0_ten = sympify_expression(expr * 10)
            self.assertEqual( (sympified_expr * 10 - prod_p0_ten).simplify() , 0)

            prod_ten_p0 = sympify_expression(10 * expr)
            self.assertEqual( (sympified_expr * 10 - prod_ten_p0).simplify() , 0)

            # other operand is sympy symbol
            prod_p0_K = sympify_expression(expr * sympify_expression('K'))
            self.assertEqual( (sympified_expr * sympify_expression('K') - prod_p0_K).simplify() , 0)

            prod_K_p0 = sympify_expression(sympify_expression('K') * expr)
            self.assertEqual( (sympified_expr * sympify_expression('K') - prod_K_p0).simplify() , 0)

    #@attr('active')
    def test_pow(self):
        pow_p0_p1 = self.p0 ** self.p1
        sympified_pow_p0_p1 = sympify_expression(pow_p0_p1)
        target_pow_p0_p1 = sympify_expression('((a*x + b*y)**exponent) ** log(c*x + d*y)')
        self.assertTrue(type(pow_p0_p1) is Pow)
        self.assertEqual( (sympified_pow_p0_p1 - target_pow_p0_p1).simplify() , 0 )

        pow_p1_p0 = self.p1 ** self.p0
        sympified_pow_p1_p0 = sympify_expression(pow_p1_p0)
        target_pow_p1_p0 = sympify_expression('log(c*x + d*y) ** ((a*x + b*y)**exponent)')
        self.assertTrue(type(pow_p1_p0) is Pow)
        self.assertEqual( (sympified_pow_p1_p0 - target_pow_p1_p0).simplify() , 0 )

        for expr, sympified_expr in zip([self.p0, self.p1], sympify_expression(['(a*x + b*y)**exponent', 'log(c*x + d*y)'])):
            # other operand is number
            pow_p0_ten = sympify_expression(expr ** 10)
            self.assertEqual( (sympified_expr ** 10 - pow_p0_ten).simplify() , 0)

            pow_ten_p0 = sympify_expression(10 ** expr)
            self.assertEqual( (10**sympified_expr - pow_ten_p0).simplify() , 0)

            # other operand is sympy symbol
            pow_p0_K = sympify_expression(expr ** sympify_expression('K'))
            self.assertEqual( (sympified_expr ** sympify_expression('K') - pow_p0_K).simplify() , 0)

            pow_K_p0 = sympify_expression(sympify_expression('K') ** expr)
            self.assertEqual( (sympify_expression('K')**sympified_expr - pow_K_p0).simplify() , 0)

class TestExpressionConverter(unittest.TestCase):
    #@attr('active')
    def test_polynomial(self):
        sympy_poly = sympify_expression('x + x**2 - a * x*y')
        poly = Expression(sympy_poly, ['x','y'])

        self.assertTrue( type(poly) is Polynomial )
        self.assertEqual( (sympify_expression(poly) - sympy_poly).simplify() , 0 )

    #@attr('active')
    def test_function(self):
        a, x = sp.symbols('a x')
        my_function = sp.Function('my_function')
        sympy_function = my_function(a*x)
        function = Expression(sympy_function, ['x'])

        self.assertEqual( (sympify_expression(function) - sympy_function).simplify() , 0 )
        self.assertTrue( type(function) is Function )
        self.assertEqual(function.symbol, 'my_function')
        self.assertEqual(len(function.arguments), 1)
        self.assertTrue( type(function.arguments[0]) is Polynomial )
        np.testing.assert_array_equal(function.arguments[0].expolist, [[1]])
        np.testing.assert_array_equal(function.arguments[0].coeffs, [a])

    #@attr('active')
    def test_follow_functions(self):
        string_expression = 'f1(x,y) ** f2(x*y,z) + f1(y,z) * f2(z) + f2(y*x,z)'
        x, y, z = sp.symbols('x y z')
        polysymbols = ['x', 'y', 'z']
        expression, functions = Expression(string_expression, polysymbols, follow_functions=True)

        expression.derive(2) # derivative by "z"

        found_f1_x_y = found_f2_xy_z = found_f1_y_z = found_f2_z = False

        self.assertEqual( (sympify_expression(expression) - sympify_expression(string_expression)).simplify() , 0 )
        self.assertEqual( len(functions) , 4 )
        for function in functions:
            self.assertTrue( type(function) is Function )
            if function.symbol == 'f1':
                if sympify_expression(function.arguments) == [x,y]:
                    self.assertEqual(function.derivative_symbols, set(['f1']))
                    found_f1_x_y = True
                elif sympify_expression(function.arguments) == [y,z]:
                    self.assertEqual(function.derivative_symbols, set(['f1','df1d1']))
                    found_f1_y_z = True
            elif function.symbol == 'f2':
                if sympify_expression(function.arguments) == [x*y,z]:
                    self.assertEqual(function.derivative_symbols, set(['f2','df2d1']))
                    found_f2_xy_z = True
                elif sympify_expression(function.arguments) == [z]:
                    self.assertEqual(function.derivative_symbols, set(['f2','df2d0']))
                    found_f2_z = True

        for found_call in [found_f1_x_y, found_f2_xy_z, found_f1_y_z, found_f2_z]:
            self.assertTrue(found_call)
