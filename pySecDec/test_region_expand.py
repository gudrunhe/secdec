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
    def test_apply_regions_small(self):
        poly = Polynomial.from_expression('t + u + u**2', ['u'])
        sector = Sector([poly])
        exp_param = 't'

        regions = list(apply_regions( exp_param, sector, workdir = 'tmpdir_test_apply_regions_python' + python_major_version ))
        print(regions)

        target_Jacobians = [sp.sympify('1'), sp.sympify(' zz')]
        target_polys = [sp.sympify('1 * (t* zz + u**2 + u)'), sp.sympify('zz**1 * (t + u**2*zz + u)')]

        try:
            for target_poly, target_Jacobian, region in zip(target_polys, target_Jacobians, regions):
                self.assertEqual( (sp.sympify(region.cast[0])-target_poly).simplify() , 0)
                self.assertEqual( (sp.sympify(region.Jacobian)-target_Jacobian).simplify() , 0)
        except AssertionError:
            for target_poly, target_Jacobian, region in zip(target_polys, target_Jacobians, reversed(regions)):
                self.assertEqual( (sp.sympify(region.cast[0])-target_poly).simplify() , 0)
                self.assertEqual( (sp.sympify(region.Jacobian)-target_Jacobian).simplify() , 0)

    # TODO: complete or remove
    #@attr('active')
    def test_apply_regions(self):
        poly_u = Polynomial.from_expression('x0 + x1', ['x0','x1'])
        poly_f = Polynomial.from_expression('+ (-psq*r)*x1**2 + (-2*psq*r - psq)*x0*x1 + (-psq*r)*x0**2', ['x0','x1'])
        sector = Sector([poly_u,poly_f])
        exp_param = 'r'

        regions = list(apply_regions( exp_param, sector, workdir = 'tmpdir_test_apply_regions_python' + python_major_version ))
