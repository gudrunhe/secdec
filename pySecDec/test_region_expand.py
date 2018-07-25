from .region_expand import *
from .decomposition.common import Sector
from .algebra import Polynomial, ExponentiatedPolynomial
from .misc import argsort_2D_array
from nose.plugins.attrib import attr
import numpy as np
import sympy as sp
import sys
import unittest

python_major_version = sys.version[0]

def sort_2D_array(array):
    'Use the .misc.argsort_2D_array function to sort an array'
    return array[argsort_2D_array(array)]

#@attr('active')
class  TestFindRegions(unittest.TestCase):
    #@attr('active')
    def test_find_regions_1L(self):
        poly_u = Polynomial.from_expression('x0 + x1', ['z','x0','x1'])
        poly_f = Polynomial.from_expression('+ (-psq*z)*x1**2 + (-2*psq*z - psq)*x0*x1 + (-psq*z)*x0**2', ['z','x0','x1'])
        exp_param_index = 0
        target_regions = sort_2D_array(np.array([[1,0,0], [1,0,-1], [1,-1,0]]))

        regions = find_regions( exp_param_index, [poly_u, poly_f], workdir = 'tmpdir_test_find_regions_1L_python' + python_major_version )

        np.testing.assert_array_equal(sort_2D_array(regions), target_regions)

    #@attr('active')
    def test_find_regions_2L(self):
        poly_u = Polynomial.from_expression(' + (1)*x0*x1 + (1)*x0*x3 + (1)*x0*x4 + (1)*x0*x5 + (1)*x1*x2 + (1)*x1*x4 + (1)*x2*x3 + (1)*x2*x4 + (1)*x2*x5 + (1)*x3*x4 + (1)*x4*x5', ['z','x0','x1','x2','x3','x4','x5'])
        poly_f = Polynomial.from_expression('+ ((z*m)**2)*x3**2*x4 + (2*(z*m)**2)*x2*x3*x4 + ((z*m)**2)*x2*x3**2 + ((z*m)**2)*x2**2*x5 + ((z*m)**2)*x2**2*x4 + ((z*m)**2)*x2**2*x3 + (s)*x1*x3*x4 + (s)*x1*x2*x4 + (s)*x1*x2*x3 + ((z*m)**2)*x1*x2**2 + (M**2)*x1**2*x4 + (M**2)*x1**2*x2 + (s)*x0*x3*x4 + ((z*m)**2)*x0*x3**2 + (s)*x0*x2*x5 + (s)*x0*x2*x4 + (s)*x0*x2*x3 + (2*M**2)*x0*x1*x4 + (s)*x0*x1*x3 + (s)*x0*x1*x2 + (M**2)*x0*x1**2 + (M**2)*x0**2*x5 + (M**2)*x0**2*x4 + (M**2)*x0**2*x3 + (M**2)*x0**2*x1', ['z','x0','x1','x2','x3','x4','x5'])
        exp_param_index = 0
        target_regions = sort_2D_array(np.array([[1,0,-2,-2,-2,-2,-4], [1,0,0,-2,-2,-2,-2], [1,0,0,-2,0,0,0],[1,0,0,0,-2,0,-2],[1,0,0,0,0,0,0]]))

        regions = find_regions( exp_param_index, [poly_u, poly_f], workdir = 'tmpdir_test_find_regions_2L_python' + python_major_version )

        np.testing.assert_array_equal(sort_2D_array(regions), target_regions)

    #@attr('active')
    def test_find_regions_small(self):
        poly = Polynomial.from_expression('z*t + u + u**2', ['z','u'])
        exp_param_index = 0
        target_regions = sort_2D_array(np.array([[1,0], [1,1]]))

        regions = find_regions( exp_param_index, [ poly], workdir = 'tmpdir_test_find_regions_small_python' + python_major_version )

        np.testing.assert_array_equal(sort_2D_array(regions), target_regions)

#@attr('active')
class TestExpansionByRegions(unittest.TestCase):
    #@attr('active')
    def test_apply_region_small(self):
        poly = Polynomial.from_expression('z*t + u + u**2', ['z','u'])
        expansion_parameter_index = 0
        regions = np.array([[1,0], [1,1]]) # [(z->z^1, u->z^0*u), (z->z^1, u->z^1*u)]

        target_polys = [sp.sympify('z*t + u + u**2'), sp.sympify('z*t + z*u + (z*u)**2')]

        for region, target_poly in zip(regions,target_polys):
            polys = apply_region([poly.copy()], region, expansion_parameter_index)
            self.assertEqual(sp.sympify(polys[0]-target_poly).simplify(), 0)

    #@attr('active')
    def test_apply_region(self):
        poly_u = Polynomial.from_expression('x0 + x1', ['z','x0','x1'])
        poly_f = Polynomial.from_expression('-s*x0*x1 + z*msq*(x0+x1)**2', ['z','x0','x1'])
        region_vectors = np.array([[1,0,0],[1,-1,0],[1,0,-1]])
        expansion_parameter_index = 0
        target_polynomial1 = [sp.sympify('x0 + x1'), sp.sympify('-s*x0*x1 + z*msq*(x0+x1)**2')]
        target_polynomial2 = [sp.sympify('z**(-1)*x0 + x1'), sp.sympify('-s*(z**(-1)*x0)*x1 + z*msq*((z**(-1)*x0)+x1)**2')]
        target_polynomial3 = [sp.sympify('x0 + (z**(-1)*x1)'), sp.sympify('-s*x0*(z**(-1)*x1) + z*msq*(x0+(z**(-1)*x1))**2')]
        target_polynomials=[target_polynomial1, target_polynomial2, target_polynomial3]

        for region_vector,target_polynomial in zip(region_vectors,target_polynomials):
            polys=apply_region([poly_u.copy(),poly_f.copy()],region_vector,expansion_parameter_index)
            print(polys)
            self.assertEqual(sp.sympify(polys[0]-target_polynomial[0]).simplify(), 0)
            self.assertEqual(sp.sympify(polys[1]-target_polynomial[1]).simplify(), 0)


    #@attr('active')
    def test_derivative_product_1(self):
        p1 = lambda expo: ExponentiatedPolynomial([[1,0,0,0],[0,2,0,0]],['a','b'],expo,['x','y','p1','p2'])
        p2 = lambda expo: ExponentiatedPolynomial([[1,2,0,0],[1,1,0,0]],['c','d'],expo,['x','y','p1','p2'])
        numerator = Polynomial([[1,2,1,1]],['e'],['x','y','p1','p2'])
        target_numerator = sp.sympify('e*y^2*p1*p2*(4*a*x*p2 + p1*(3*x*y*(d + c*y) + p2))')
        derive_prod_output = sp.sympify(derive_prod([p1(3),p2(2)],numerator,0,[2,3]))

        self.assertEqual(len(derive_prod_output[0]), 2)
        self.assertEqual(sp.sympify(derive_prod_output[0][0]-p1(2)).simplify(), 0)
        self.assertEqual(sp.sympify(derive_prod_output[0][1]-p2(1)).simplify(), 0)
        self.assertEqual(sp.sympify(derive_prod_output[1]-target_numerator).simplify(), 0)

    #@attr('active')
    def test_derivative_product_2(self):
        p1 = lambda expo: ExponentiatedPolynomial([[1,2,0,0,0]],['a'],expo,['x','y','p1','p2','p3'])
        p2 = lambda expo: ExponentiatedPolynomial([[3,2,0,0,0],[1,1,0,0,0]],['b','c'],expo,['x','y','p1','p2','p3'])
        p3 = lambda expo: ExponentiatedPolynomial([[2,2,0,0,0],[1,1,0,0,0]],['d','e'],expo,['x','y','p1','p2','p3'])
        numerator = Polynomial([[1,2,2,1,1]],['f'],['x','y','p1','p2','p3'])
        target_numerator = sp.sympify('f*x*y*p1^2*p2*p3*(10*a*x*y^2*p2*p3 + p1*(3*x*y*(c + 2*b*x^2*y)*p3 + p2*(3*x*y*(e + 2*d*x*y) +2*p3)))')

        derive_prod_output = sp.sympify(derive_prod([p1(3),p2(2),p3(2)],numerator,1,[2,3,4]))

        self.assertEqual(len(derive_prod_output[0]), 3)
        self.assertEqual(sp.sympify(derive_prod_output[0][0]-p1(2)).simplify(), 0)
        self.assertEqual(sp.sympify(derive_prod_output[0][1]-p2(1)).simplify(), 0)
        self.assertEqual(sp.sympify(derive_prod_output[0][2]-p3(1)).simplify(), 0)
        self.assertEqual(sp.sympify(derive_prod_output[1]-target_numerator).simplify(), 0)

    #@attr('active')
    def test_derivative_product_3(self):
        p1 = lambda expo: ExponentiatedPolynomial([[1,0,0,0],[0,2,0,0]],['a','b'],expo,['x','y','p1','p2'])
        p2 = lambda expo: ExponentiatedPolynomial([[1,2,0,0],[1,1,0,0]],['c','d'],expo,['x','y','p1','p2'])
        numerator = Polynomial([[1,2,1,1]],['e'],['x','y','p1','p2'])
        target_numerator = sp.sympify('e*y^2*p1*p2*(2*a*x*p2 + p1*(2*x*y*(d + c*y) + p2))')

        derive_prod_output = sp.sympify(derive_prod([p1(1),p2(1)],numerator,0,[2,3]))

        self.assertEqual(len(derive_prod_output[0]), 2)
        self.assertEqual(sp.sympify(derive_prod_output[0][0]-p1(0)).simplify(), 0)
        self.assertEqual(sp.sympify(derive_prod_output[0][1]-p2(0)).simplify(), 0)
        self.assertEqual(sp.sympify(derive_prod_output[1]-target_numerator).simplify(), 0)

    #@attr('active')
    def test_derivative_product_exponent_eps(self):
        p1 = lambda expo: ExponentiatedPolynomial([[1,0,0,0],[0,2,0,0]],['a','b'],expo,['x','y','p1','p2'])
        p2 = lambda expo: ExponentiatedPolynomial([[1,2,0,0],[1,1,0,0]],['c','d'],expo,['x','y','p1','p2'])
        numerator = Polynomial([[1,2,0,0]],['e'],['x','y','p1','p2'])
        target_numerator = sp.sympify('e*y^2*(-(a*(-3 + eps)*x*p2) + p1*(-((-4 + eps)*x*y*(d + c*y)) +p2))')

        poly_list, numerator = sp.sympify(derive_prod([p1('3-eps'),p2('4-eps')],numerator,0,[2,3]))

        self.assertEqual(sp.sympify(poly_list[0]-p1('2-eps')).simplify(), 0)
        self.assertEqual(sp.sympify(poly_list[1]-p2('3-eps')).simplify(), 0)
        self.assertEqual(sp.sympify(numerator-target_numerator).simplify(), 0)
