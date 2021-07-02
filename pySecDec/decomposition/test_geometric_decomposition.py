from .geometric import *
from .common import Sector
from ..algebra import Polynomial, ExponentiatedPolynomial, LogOfPolynomial
from ..misc import argsort_2D_array, sympify_expression
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
class TestGeomethod(unittest.TestCase):
    def setUp(self):
        self.p0 = Polynomial.from_expression('x0+x1+x0*x1', ['x0','x1'])
        self.p1 = Polynomial.from_expression('1+x0+x1', ['x0','x1'])
        self.target_hull = np.array([[2,1],
                                     [1,2],
                                     [2,0],
                                     [1,0],
                                     [0,2],
                                     [0,1]])
        self.sorted_target_hull = sort_2D_array(self.target_hull)
        self.target_fan_p01 = [[(1,-1),(1,0),(0,1)],
                               [(-1,1),(1,0),(0,1)]]
        self.sorted_target_fan_p01 = sorted([sorted(cone) for cone in self.target_fan_p01])
        self.p2 = Polynomial.from_expression('x0+x0*x2+x1*x2', ['x0','x1','x2'])
        self.target_fan_p2 = [[(1,-1,-1),(1,-1,0),(1,0,0),(0,1,0),(0,0,1)],
                              [(-1,1,1),(1,0,0),(0,1,0),(0,0,1)]]
        self.sorted_target_fan_p2 = sorted([sorted(cone) for cone in self.target_fan_p2])

    #@attr('active')
    def test_Cheng_Wu(self):
        Feynman_parameters = ['x1', 'x2', 'x3']
        x1, x2, x3 = sp.symbols(Feynman_parameters)
        poly_sympy = sympify_expression('x1*x2 + 5*x2*x3')
        other_sympy = sympify_expression('a + x3*6')
        poly = Polynomial.from_expression(poly_sympy, Feynman_parameters)
        other = Polynomial.from_expression(other_sympy, Feynman_parameters)
        sector = Sector([poly], [other])

        primary_x3 = Cheng_Wu(sector)
        # ``x2`` factorizes in ``poly``
        self.assertEqual( (sympify_expression(primary_x3.cast[0].factors[0]) - x2).simplify() , 0 )
        self.assertEqual( (sympify_expression(primary_x3.cast[0].factors[1]) - (poly_sympy/x2).subs('x3',1)).simplify() , 0 )
        self.assertEqual( (sympify_expression(primary_x3.other[0]) - other_sympy.subs('x3',1)).simplify() , 0 )

        primary_x1 = Cheng_Wu(sector,0)
        # ``x2`` factorizes in ``poly``
        self.assertEqual( (sympify_expression(primary_x1.cast[0].factors[0]) - x2).simplify() , 0 )
        self.assertEqual( (sympify_expression(primary_x1.cast[0].factors[1]) - (poly_sympy/x2).subs('x1',1)).simplify() , 0 )
        self.assertEqual( (sympify_expression(primary_x1.other[0]) - other_sympy.subs('x1',1)).simplify() , 0 )

    #@attr('active')
    def test_Cheng_Wu_one_variable(self):
        U = Polynomial([[1]], [  1  ], ['x0'])
        F = Polynomial([[2]], ['msq'], ['x0'])
        initial_sector = Sector([U,F])
        primary_sector = Cheng_Wu(initial_sector)

        target_decomposed_U = sympify_expression(1)
        target_decomposed_F = sympify_expression('msq')

        self.assertEqual(  ( sympify_expression(primary_sector.cast[0]) - target_decomposed_U ).simplify() , 0   )
        self.assertEqual(  ( sympify_expression(primary_sector.cast[1]) - target_decomposed_F ).simplify() , 0   )

    def test_convex_hull(self):
        hull = convex_hull(self.p0, self.p1)

        # The ordering is not important but must be fixed to compare the arrays
        sorted_hull = sort_2D_array(hull)
        np.testing.assert_array_equal(sorted_hull, self.sorted_target_hull)

    def test_generate_fan(self):
        fan_p01 = generate_fan(self.p0,self.p1)
        for cone, target_cone in zip(fan_p01, self.sorted_target_fan_p01):
            np.testing.assert_array_equal(cone, target_cone)

        fan_p2 = generate_fan(self.p2)
        for cone, target_cone in zip(fan_p2, self.sorted_target_fan_p2):
            np.testing.assert_array_equal(cone, target_cone)

        fan_p = generate_fan(Polynomial.from_expression('1+x', ['x']))
        target_fan_p = [[[1]]]
        self.assertEqual(len(fan_p), len(target_fan_p))
        for cone, target_cone in zip(fan_p, target_fan_p):
            np.testing.assert_array_equal(cone, target_cone)

    #@attr('active')
    def test_convex_hull_exponentiated_polynomial(self):
        p0 = ExponentiatedPolynomial(self.p0.expolist, self.p0.coeffs, polysymbols=self.p0.polysymbols, exponent='8-3*eps')
        p1 = ExponentiatedPolynomial(self.p1.expolist, self.p1.coeffs, polysymbols=self.p1.polysymbols, exponent='1-2*eps')

        hull = convex_hull(p0, p1)

        # The ordering is not important but must be fixed to compare the arrays
        sorted_hull = sort_2D_array(hull)
        np.testing.assert_array_equal(sorted_hull, self.sorted_target_hull)

    #@attr('active')
    def test_triangulate(self):
        # basic consistency checks working?
        simplicial_cone = [[ 1,  0,  0], [ 0,  1,  0], [ 0, -1, -1]]
        self.assertRaisesRegexp(ValueError, 'simplicial.*already', triangulate, simplicial_cone)
        wrong_dimensionality = [ 1,  0,  0]
        self.assertRaisesRegexp(AssertionError, '(M|m)ust.*two.*dim', triangulate, wrong_dimensionality)
        two_rays = [[ 1,  0,  0], [ 0,  1,  0]]
        self.assertRaisesRegexp(AssertionError, '(M|m)ust.*at least.*dim', triangulate, two_rays)


        cone = [[ 1,  0,  0], [ 0,  1,  0], [ 0, -1, -1], [-1,  0, -1]]
        cone_normal = [[ -1, 1, 1], [ 1, 0, 0], [ 0, 1, 0], [ 0, 0, 1]]

        # useful error message?
        self.assertRaisesRegexp(
                                    OSError, 'No such file or directory.*nonexistentNormalizExecutable',
                                    triangulate, cone, normaliz='nonexistentNormalizExecutable',
                                    workdir='tmpdir_test_triangulate_python' + python_major_version
                               )

        triangulated_cones = triangulate(cone, workdir='tmpdir_test_triangulate_python' + python_major_version)
        triangulated_cones_normal = triangulate(cone_normal, workdir='tmpdir_test_triangulate_python' + python_major_version, switch_representation=True)

        # there are two possibilities for the triangualtion
        target_triangulated_cones1 = np.array([
                                                [[ 1,  0,  0], [ 0,  1,  0], [-1,  0, -1]],
                                                [[ 1,  0,  0], [ 0, -1, -1], [-1,  0, -1]]
                                            ])
        target_triangulated_cones2 = np.array([
                                                [[ 0, -1, -1], [ 0,  1,  0], [ 1,  0,  0]],
                                                [[ 0, -1, -1], [ 0,  1,  0], [-1,  0, -1]]
                                            ])
        target_triangulated_cones1_normal = np.array([
                                                [[ 1,  1,  0], [ 1,  0,  1], [ 0,  0, 1]],
                                                [[ 1,  1,  0], [ 0, 1, 0], [0,  0, 1]]
                                            ])
        target_triangulated_cones2_normal = np.array([
                                                [[ 1, 0, 1], [ 0,  1,  0], [ 0,  0,  1]],
                                                [[ 1, 0, 1], [ 0,  1,  0], [ 1,  1, 0]]
                                            ])

        # should get one of these triangulations
        # The ordering is not important but must be fixed to compare the arrays
        try:
            np.testing.assert_array_equal(sort_2D_array(triangulated_cones[0]), sort_2D_array(target_triangulated_cones1[0]))
            np.testing.assert_array_equal(sort_2D_array(triangulated_cones[1]), sort_2D_array(target_triangulated_cones1[1]))
        except AssertionError:
            try:
                np.testing.assert_array_equal(sort_2D_array(triangulated_cones[0]), sort_2D_array(target_triangulated_cones1[1]))
                np.testing.assert_array_equal(sort_2D_array(triangulated_cones[1]), sort_2D_array(target_triangulated_cones1[0]))
            except AssertionError:
                try:
                    np.testing.assert_array_equal(sort_2D_array(triangulated_cones[0]), sort_2D_array(target_triangulated_cones2[0]))
                    np.testing.assert_array_equal(sort_2D_array(triangulated_cones[1]), sort_2D_array(target_triangulated_cones2[1]))
                except:
                    np.testing.assert_array_equal(sort_2D_array(triangulated_cones[0]), sort_2D_array(target_triangulated_cones2[1]))
                    np.testing.assert_array_equal(sort_2D_array(triangulated_cones[1]), sort_2D_array(target_triangulated_cones2[0]))

        try:
            np.testing.assert_array_equal(sort_2D_array(triangulated_cones_normal[0]), sort_2D_array(target_triangulated_cones1_normal[0]))
            np.testing.assert_array_equal(sort_2D_array(triangulated_cones_normal[1]), sort_2D_array(target_triangulated_cones1_normal[1]))
        except AssertionError:
            try:
                np.testing.assert_array_equal(sort_2D_array(triangulated_cones_normal[0]), sort_2D_array(target_triangulated_cones1_normal[1]))
                np.testing.assert_array_equal(sort_2D_array(triangulated_cones_normal[1]), sort_2D_array(target_triangulated_cones1_normal[0]))
            except AssertionError:
                try:
                    np.testing.assert_array_equal(sort_2D_array(triangulated_cones_normal[0]), sort_2D_array(target_triangulated_cones2_normal[0]))
                    np.testing.assert_array_equal(sort_2D_array(triangulated_cones_normal[1]), sort_2D_array(target_triangulated_cones2_normal[1]))
                except:
                    np.testing.assert_array_equal(sort_2D_array(triangulated_cones_normal[0]), sort_2D_array(target_triangulated_cones2_normal[1]))
                    np.testing.assert_array_equal(sort_2D_array(triangulated_cones_normal[1]), sort_2D_array(target_triangulated_cones2_normal[0]))

    def test_triangulate_1D(self):
        cone = [[1]]
        target_cone = np.array([[[1]]])
        triangulated_cone = triangulate(cone, workdir='tmpdir_test_triangulate_1D_python' + python_major_version, switch_representation=True)
        np.testing.assert_array_equal(triangulated_cone, target_cone)

    #@attr('active')
    def test_transform_variables(self):
        x0 = Polynomial.from_expression('x0',['x0','x1','x2'])
        x1 = Polynomial.from_expression('x1',['x0','x1','x2'])
        x2 = Polynomial.from_expression('x2',['x0','x1','x2'])

        y0 = Polynomial.from_expression('y0',['y0','y1','y2','y3'])
        y1 = Polynomial.from_expression('y1',['y0','y1','y2','y3'])
        y2 = Polynomial.from_expression('y2',['y0','y1','y2','y3'])
        y3 = Polynomial.from_expression('y3',['y0','y1','y2','y3'])

        composite_polynomial = x0 * x0 * x1
        exponentiated_polynomial = ExponentiatedPolynomial([[2,1,0],[0,0,0]], [1,2], polysymbols='x', exponent='exponent')
        log_of_polynomial = LogOfPolynomial([[2,1,0],[0,1,1]], [1,2], polysymbols='x')

        transformation = np.array([[ 0, 0, 0, 1],   # x0 -> y3
                                   [ 1, 1, 0, 0],   # x1 -> y0*y1
                                   [ 1, 1, 1,-2]])  # x2 -> y0*y1*y2*y3**-2

        transformed_x0 = y3
        transformed_x1 = y0 * y1
        transformed_x2 = Polynomial(expolist=[[1, 1, 1,-2]], coeffs=[1], polysymbols='y')
        transformed_composite_polynomial = transformed_x0 * transformed_x0 * transformed_x1
        transformed_exponentiated_polynomial = ExponentiatedPolynomial([[1,1,0,2],[0,0,0,0]], [1,2], polysymbols='y', exponent='exponent')
        transformed_log_of_polynomial = LogOfPolynomial([[1,1,0,2],[2,2,1,-2]], [1,2], polysymbols='y')

        self.assertEqual( ( sympify_expression(transformed_x0) - sympify_expression(transform_variables(x0, transformation)) ).simplify() , 0)
        self.assertEqual( (sympify_expression(transformed_x1) - sympify_expression(transform_variables(x1, transformation))).simplify() , 0)
        self.assertEqual( (sympify_expression(transformed_x2) - sympify_expression(transform_variables(x2, transformation))).simplify() , 0)
        self.assertEqual( (sympify_expression(transformed_composite_polynomial) - sympify_expression(transform_variables(composite_polynomial, transformation))).simplify() , 0)

        self.assertTrue(type(transform_variables(exponentiated_polynomial, transformation)) is ExponentiatedPolynomial)
        self.assertEqual( (sympify_expression(transformed_exponentiated_polynomial) - sympify_expression(transform_variables(exponentiated_polynomial, transformation))).simplify() , 0)

        self.assertTrue(type(transform_variables(log_of_polynomial, transformation)) is LogOfPolynomial)
        self.assertEqual( (sympify_expression(transformed_log_of_polynomial) - sympify_expression(transform_variables(log_of_polynomial, transformation))).simplify() , 0)

    #@attr('active')
    def test_2D_geometric_decomposition(self):
        poly = Polynomial.from_expression('x1 + x2 + x1*x2', ['dummy','x1','x2'])
        sector = Sector([poly])
        indices = [1,2]
        subsectors = list( geometric_decomposition(sector, indices, workdir='tmpdir_test_2D_geometric_decomposition_python' + python_major_version) )

        target_general_Jacobian = sympify_expression('x1**-2 * x2**-2 * x3')
        target_general_poly = sympify_expression('x1**-1 * x2**-1 * x3 * (x1 + x2 + x3)')

        self.assertEqual(len(subsectors), 3)

        for i,subsector in enumerate(subsectors):
            Jacobian = sympify_expression(subsector.Jacobian)
            poly = sympify_expression(subsector.cast[0])

            target_Jacobian = target_general_Jacobian.subs('x%i'%(i+1), 1)
            target_poly = target_general_poly.subs('x%i'%(i+1), 1)
            for j in range(i+1,3+1):
                target_Jacobian = target_Jacobian.subs('x%i'%j, 'x%i'%(j-1))
                target_poly = target_poly.subs('x%i'%j, 'x%i'%(j-1))

            self.assertEqual( (poly-target_poly).simplify() , 0)
            self.assertEqual( (Jacobian-target_Jacobian).simplify() , 0)

    #@attr('active')
    def test_2D_geometric_decomposition_ku(self):
        poly = Polynomial.from_expression('A*x1 + B*x2 + C*x1*x2', ['dummy','x1','x2'])
        sector = Sector([poly])
        indices = [1,2]
        subsectors = list( geometric_decomposition_ku(sector, indices, workdir='tmpdir_test_2D_geometric_decomposition_ku_python' + python_major_version) )
        print(subsectors)
        self.assertEqual(len(subsectors), 2)

        target_Jacobians = [sympify_expression('x2**1 '), sympify_expression('x1**1 ')]
        target_polys = [sympify_expression('x2**1 * (A*x1 + B + C*x1*x2)'), sympify_expression('x1**1 * (A + B*x2 + C*x1*x2)')]

        try:
            for target_poly, target_Jacobian, subsector in zip(target_polys, target_Jacobians,subsectors):
                try:
                    self.assertEqual( (sympify_expression(subsector.cast[0])-target_poly).simplify() , 0)
                    self.assertEqual( (sympify_expression(subsector.Jacobian)-target_Jacobian).simplify() , 0)
                except AssertionError:
                    self.assertEqual( (sympify_expression(subsector.cast[0])-target_poly.subs([('x1','x2'),('x2','x1')],simultaneous=True)).simplify() , 0)
                    self.assertEqual( (sympify_expression(subsector.Jacobian)-target_Jacobian.subs([('x1','x2'),('x2','x1')],simultaneous=True)).simplify() , 0)
        except AssertionError:
            for target_poly, target_Jacobian, subsector in zip(target_polys, target_Jacobians,reversed(subsectors)):
                try:
                    self.assertEqual( (sympify_expression(subsector.cast[0])-target_poly).simplify() , 0)
                    self.assertEqual( (sympify_expression(subsector.Jacobian)-target_Jacobian).simplify() , 0)
                except AssertionError:
                    self.assertEqual( (sympify_expression(subsector.cast[0])-target_poly.subs([('x1','x2'),('x2','x1')],simultaneous=True)).simplify() , 0)
                    self.assertEqual( (sympify_expression(subsector.Jacobian)-target_Jacobian.subs([('x1','x2'),('x2','x1')],simultaneous=True)).simplify() , 0)

    #@attr('active')
    def test_geometric_ku_lower_dimensional_cones(self):
        polysymbols = ['x0','x1','x3','x4','x5']
        poly1 = Polynomial.from_expression(
                                               '''
                                                   + (4) + (2)*x5 + (9/5)*x4 + (-7/10)*x3 + (-7/40)*x3*x5 + (-63/400)*x3*x4 + (-4/5)*x1
                                                   + (-1/5)*x1*x5 + (-9/50)*x1*x4 + (7/50)*x1*x3 + (-7/5)*x0 + (-7/20)*x0*x5
                                                   + (-63/200)*x0*x4 + (49/200)*x0*x3
                                               ''',
                                               polysymbols
                                          )
        poly2 = Polynomial.from_expression(
                                              '''
                                                   + (2) + (1)*x5 + (9/10)*x4 + (-7/4)*x3 + (-7/10)*x3*x5 + (-63/80)*x3*x4
                                                   + (49/200)*x3**2 + (49/800)*x3**2*x5 + (441/8000)*x3**2*x4 + (9/50)*x1*x4
                                                   + (7/25)*x1*x3 + (7/100)*x1*x3*x5 + (63/1000)*x1*x3*x4 + (-49/1000)*x1*x3**2
                                                   + (7/20)*x0*x5 + (49/100)*x0*x3 + (441/4000)*x0*x3*x4 + (-343/4000)*x0*x3**2
                                                   + (-14/25)*x0*x1 + (-7/50)*x0*x1*x5 + (-63/500)*x0*x1*x4 + (49/500)*x0*x1*x3
                                              ''',
                                              polysymbols
                                         )
        sector = Sector([poly1,poly2])

        # should not error
        subsectors = list( geometric_decomposition_ku(sector, workdir='tmpdir_test_geometric_ku_lower_dimensional_cones_python' + python_major_version) )

    #@attr('active')
    def test_3D_geometric_decomposition(self):
        # 3D test case where triangulation is needed
        poly = Polynomial.from_expression('A*1 + B*x1 + C*x2 + D*x3 + E*x1*x2', ['x1','x2','x3']) # pyramid
        sector = Sector([poly])
        subsectors = list( geometric_decomposition(sector, workdir='tmpdir_test_3D_geometric_decomposition_python' + python_major_version) )

    #@attr('active')
    def test_3D_geometric_decomposition_selected_indices(self):
        # 3D test case where triangulation is needed
        poly = Polynomial.from_expression('A*1 + B*x1 + C*x2 + D*x3 + E*x1*x2', ['x1','dummy','x2','x3']) # pyramid
        sector = Sector([poly])
        indices = [0,2,3]
        subsectors = list( geometric_decomposition(sector, indices, workdir='tmpdir_test_3D_geometric_decomposition_selected_indices_python' + python_major_version) )
