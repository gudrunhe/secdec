from .polytope import *
from .misc import argsort_2D_array
from nose.plugins.attrib import attr
import sys
import unittest

python_major_version = sys.version[0]

def sort_2D_array(array):
    'Use the .misc.argsort_2D_array function to sort an array'
    return array[argsort_2D_array(array)]

#@attr('active')
class TestPolytope(unittest.TestCase):
    def setUp(self):
        self.vertices = [[2,1],
                         [1,2],
                         [2,0],
                         [1,0],
                         [0,2],
                         [0,1]]

        self.vertices_with_inside = [[1,1]] + self.vertices

        self.facets = [[ 0, 1, 0],
                       [ 1, 0, 0],
                       [-1, 0, 2],
                       [ 0,-1, 2],
                       [ 1, 1,-1],
                       [-1,-1, 3]]

        self.facets_with_outside = [[ 0, 1, 5]] + self.facets

    def test_init(self):
        # must provide either facets or vertices
        self.assertRaisesRegexp(TypeError, 'either.*vertices.*or.*facets', Polytope, facets=self.facets, vertices=self.vertices)
        self.assertRaisesRegexp(TypeError, 'either.*vertices.*or.*facets', Polytope)

        inarray = np.array([[1,2,3],
                            [4,5,6]])

        polytope1 = Polytope(facets=inarray)
        polytope2 = Polytope(vertices=inarray)

        inarray[0,0] = 10

        self.assertEqual(polytope1.facets[0,0], 1)
        self.assertEqual(polytope2.vertices[0,0], 1)

    #@attr('active')
    @attr('slow')
    def test_large_array(self):
        # numpy abbreviates long output by default, but we do not want that
        from numpy import array
        long_array = \
        array([[0, 0, 0, 2, 1, 1],
               [0, 0, 0, 2, 1, 3],
               [0, 0, 0, 2, 3, 1],
               [0, 0, 0, 3, 1, 3],
               [0, 0, 0, 3, 3, 1],
               [0, 0, 0, 4, 1, 1],
               [0, 0, 0, 4, 1, 3],
               [0, 0, 0, 4, 2, 3],
               [0, 0, 0, 4, 3, 1],
               [0, 0, 0, 4, 3, 2],
               [0, 0, 1, 2, 1, 0],
               [0, 0, 1, 2, 3, 0],
               [0, 0, 1, 3, 3, 0],
               [0, 0, 1, 3, 4, 1],
               [0, 0, 1, 4, 1, 0],
               [0, 0, 1, 4, 3, 0],
               [0, 0, 2, 0, 1, 1],
               [0, 0, 2, 0, 1, 3],
               [0, 0, 2, 0, 3, 1],
               [0, 0, 2, 1, 1, 0],
               [0, 0, 2, 1, 3, 0],
               [0, 0, 2, 3, 4, 0],
               [0, 0, 2, 4, 3, 0],
               [0, 0, 3, 0, 1, 1],
               [0, 0, 3, 0, 1, 3],
               [0, 0, 3, 0, 3, 1],
               [0, 0, 3, 1, 1, 0],
               [0, 0, 3, 1, 1, 3],
               [0, 0, 3, 1, 2, 3],
               [0, 0, 3, 1, 3, 0],
               [0, 0, 3, 1, 4, 1],
               [0, 0, 3, 2, 4, 0],
               [0, 0, 3, 4, 1, 0],
               [0, 0, 3, 4, 2, 0],
               [0, 1, 0, 2, 0, 1],
               [0, 1, 0, 2, 0, 3],
               [0, 1, 0, 3, 0, 3],
               [0, 1, 0, 3, 1, 4],
               [0, 1, 0, 4, 0, 1],
               [0, 1, 0, 4, 0, 3],
               [0, 1, 1, 1, 0, 0],
               [0, 1, 1, 4, 0, 0],
               [0, 1, 2, 0, 0, 0],
               [0, 1, 2, 0, 0, 3],
               [0, 1, 2, 0, 3, 0],
               [0, 1, 3, 0, 0, 0],
               [0, 1, 3, 0, 0, 3],
               [0, 1, 3, 0, 1, 4],
               [0, 1, 3, 0, 3, 0],
               [0, 1, 3, 0, 4, 1],
               [0, 1, 3, 1, 0, 3],
               [0, 1, 3, 4, 0, 0],
               [0, 2, 0, 0, 1, 1],
               [0, 2, 0, 0, 1, 3],
               [0, 2, 0, 0, 3, 1],
               [0, 2, 0, 1, 0, 1],
               [0, 2, 0, 1, 0, 3],
               [0, 2, 0, 3, 0, 4],
               [0, 2, 0, 4, 0, 3],
               [0, 2, 1, 0, 0, 0],
               [0, 2, 1, 0, 0, 3],
               [0, 2, 1, 0, 3, 0],
               [0, 2, 3, 0, 0, 0],
               [0, 2, 3, 0, 0, 4],
               [0, 2, 3, 0, 4, 0],
               [0, 2, 3, 4, 0, 0],
               [0, 3, 0, 0, 1, 1],
               [0, 3, 0, 0, 1, 3],
               [0, 3, 0, 0, 3, 1],
               [0, 3, 0, 1, 0, 1],
               [0, 3, 0, 1, 0, 3],
               [0, 3, 0, 1, 1, 4],
               [0, 3, 0, 1, 3, 1],
               [0, 3, 0, 1, 3, 2],
               [0, 3, 0, 2, 0, 4],
               [0, 3, 0, 4, 0, 1],
               [0, 3, 0, 4, 0, 2],
               [0, 3, 1, 0, 0, 0],
               [0, 3, 1, 0, 0, 3],
               [0, 3, 1, 0, 1, 4],
               [0, 3, 1, 0, 3, 0],
               [0, 3, 1, 0, 4, 1],
               [0, 3, 1, 1, 3, 0],
               [0, 3, 1, 4, 0, 0],
               [0, 3, 2, 0, 0, 0],
               [0, 3, 2, 0, 0, 4],
               [0, 3, 2, 0, 4, 0],
               [0, 3, 2, 4, 0, 0],
               [1, 0, 0, 2, 0, 0],
               [1, 0, 0, 2, 0, 2],
               [1, 0, 0, 2, 2, 0],
               [1, 0, 0, 3, 0, 2],
               [1, 0, 0, 3, 2, 0],
               [1, 0, 0, 4, 0, 0],
               [1, 0, 0, 4, 0, 2],
               [1, 0, 0, 4, 2, 0],
               [1, 0, 0, 4, 2, 2],
               [1, 0, 1, 1, 0, 0],
               [1, 0, 1, 1, 0, 2],
               [1, 0, 1, 1, 3, 0],
               [1, 0, 1, 2, 3, 0],
               [1, 0, 2, 1, 0, 0],
               [1, 0, 2, 1, 0, 2],
               [1, 0, 2, 1, 3, 0],
               [1, 0, 2, 2, 0, 2],
               [1, 0, 2, 2, 4, 0],
               [1, 0, 2, 4, 0, 0],
               [1, 0, 2, 4, 2, 0],
               [1, 1, 0, 1, 0, 0],
               [1, 1, 0, 1, 0, 3],
               [1, 1, 0, 1, 2, 0],
               [1, 1, 0, 2, 0, 3],
               [1, 1, 1, 0, 0, 0],
               [1, 1, 1, 0, 0, 3],
               [1, 1, 1, 0, 3, 0],
               [1, 1, 2, 0, 0, 0],
               [1, 1, 2, 0, 0, 3],
               [1, 1, 2, 0, 3, 0],
               [1, 2, 0, 1, 0, 0],
               [1, 2, 0, 1, 0, 3],
               [1, 2, 0, 1, 2, 0],
               [1, 2, 0, 2, 0, 4],
               [1, 2, 0, 2, 2, 0],
               [1, 2, 0, 4, 0, 0],
               [1, 2, 0, 4, 0, 2],
               [1, 2, 1, 0, 0, 0],
               [1, 2, 1, 0, 0, 3],
               [1, 2, 1, 0, 3, 0],
               [1, 2, 2, 0, 0, 0],
               [1, 2, 2, 0, 0, 4],
               [1, 2, 2, 0, 4, 0],
               [1, 2, 2, 4, 0, 0],
               [2, 0, 0, 0, 1, 1],
               [2, 0, 0, 0, 1, 3],
               [2, 0, 0, 0, 3, 1],
               [2, 0, 0, 1, 0, 0],
               [2, 0, 0, 1, 0, 2],
               [2, 0, 0, 1, 1, 3],
               [2, 0, 0, 1, 2, 0],
               [2, 0, 0, 1, 3, 1],
               [2, 0, 0, 2, 0, 2],
               [2, 0, 0, 2, 1, 3],
               [2, 0, 0, 2, 2, 0],
               [2, 0, 0, 2, 2, 3],
               [2, 0, 0, 2, 3, 1],
               [2, 0, 0, 2, 3, 2],
               [2, 0, 0, 3, 0, 2],
               [2, 0, 0, 3, 2, 0],
               [2, 0, 0, 3, 2, 2],
               [2, 0, 0, 4, 0, 0],
               [2, 0, 0, 4, 0, 1],
               [2, 0, 0, 4, 1, 0],
               [2, 0, 0, 4, 1, 1],
               [2, 0, 1, 0, 1, 1],
               [2, 0, 1, 0, 1, 3],
               [2, 0, 1, 0, 3, 1],
               [2, 0, 1, 1, 0, 0],
               [2, 0, 1, 1, 0, 2],
               [2, 0, 1, 1, 1, 3],
               [2, 0, 1, 1, 2, 0],
               [2, 0, 1, 1, 2, 3],
               [2, 0, 1, 1, 4, 1],
               [2, 0, 1, 2, 0, 2],
               [2, 0, 1, 2, 3, 0],
               [2, 0, 1, 3, 2, 1],
               [2, 0, 1, 4, 0, 0],
               [2, 0, 1, 4, 1, 0],
               [2, 1, 0, 0, 1, 1],
               [2, 1, 0, 0, 1, 3],
               [2, 1, 0, 0, 3, 1],
               [2, 1, 0, 1, 0, 0],
               [2, 1, 0, 1, 0, 2],
               [2, 1, 0, 1, 1, 4],
               [2, 1, 0, 1, 2, 0],
               [2, 1, 0, 1, 3, 1],
               [2, 1, 0, 1, 3, 2],
               [2, 1, 0, 2, 0, 3],
               [2, 1, 0, 2, 2, 0],
               [2, 1, 0, 3, 1, 2],
               [2, 1, 0, 4, 0, 0],
               [2, 1, 0, 4, 0, 1],
               [2, 1, 1, 0, 1, 1],
               [2, 1, 1, 0, 1, 4],
               [2, 1, 1, 0, 4, 1],
               [2, 1, 1, 1, 0, 0],
               [2, 1, 1, 1, 0, 3],
               [2, 1, 1, 1, 3, 0],
               [2, 1, 1, 3, 1, 1],
               [2, 1, 1, 4, 0, 0]])

        # error if output abbreviated by numpy
        Polytope(vertices=long_array).complete_representation(workdir='tmpdir_test_large_array_vertices_python' + python_major_version)
        Polytope(facets=long_array).complete_representation(workdir='tmpdir_test_large_array_facets_python' + python_major_version)

    #@attr('active')
    def test_vertex_incidence_lists(self):
        polytopes = [
                        Polytope(vertices=self.vertices_with_inside),
                        Polytope(facets=self.facets_with_outside),
                        Polytope(facets=self.facets),
                        Polytope(vertices=self.vertices)
                    ]

        incidences = []
        for polytope in polytopes:
            # sensible error message if `complete_representaion` not run before?
            self.assertRaisesRegexp(AssertionError, 'complete_representation.*first', polytope.vertex_incidence_lists)

            polytope.complete_representation(workdir='tmpdir_test_vertex_incidence_lists_python' + python_major_version)
            incidences.append(polytope.vertex_incidence_lists())

        target_incidence_lists = {
                                     (1,0): np.array([[ 0, 1, 0], [ 1, 1,-1]]),
                                     (0,1): np.array([[ 1, 0, 0], [ 1, 1,-1]]),
                                     (2,0): np.array([[ 0, 1, 0], [-1, 0, 2]]),
                                     (0,2): np.array([[ 1, 0, 0], [ 0,-1, 2]]),
                                     (2,1): np.array([[-1, 0, 2], [-1,-1, 3]]),
                                     (1,2): np.array([[ 0,-1, 2], [-1,-1, 3]])
                                 }
        for polytope, incidence in zip(polytopes,incidences):
            self.assertEqual(len(incidence.keys()), 6)
            for key, value in incidence.items():
                self.assertTrue(key in target_incidence_lists.keys())

                # The ordering is not important but must be fixed to compare the arrays
                np.testing.assert_array_equal(sort_2D_array(polytope.facets[value]), sort_2D_array(target_incidence_lists[key]))

    #@attr('active')
    def test_vertex2facet(self):
        polytope1 = Polytope(vertices=self.vertices)
        polytope2 = Polytope(vertices=self.vertices_with_inside)

        # useful error message?
        for polytope in (polytope1, polytope2):
            self.assertRaisesRegexp(
                                        OSError, 'No such file or directory.*nonexistentNormalizExecutable',
                                        polytope.complete_representation, normaliz='nonexistentNormalizExecutable',
                                        workdir='tmpdir_test_vertex2facet_python' + python_major_version
                                   )

        polytope1.complete_representation(workdir='tmpdir1_test_vertex2facet_python' + python_major_version)
        polytope2.complete_representation(workdir='tmpdir2_test_vertex2facet_python' + python_major_version)

        self.assertRaisesRegexp(ValueError, '(B|b)oth.*already', polytope1.complete_representation, workdir='tmpdir3_test_vertex2facet_python' + python_major_version)
        self.assertRaisesRegexp(ValueError, '(B|b)oth.*already', polytope2.complete_representation, workdir='tmpdir4_test_vertex2facet_python' + python_major_version)

        # The ordering is not important but must be fixed to compare the arrays
        np.testing.assert_array_equal( sort_2D_array(np.array(polytope1.vertices)), sort_2D_array(np.array(self.vertices)) )
        np.testing.assert_array_equal( sort_2D_array(np.array(polytope2.vertices)), sort_2D_array(np.array(self.vertices)) )
        np.testing.assert_array_equal( sort_2D_array(np.array(polytope1.facets)), sort_2D_array(np.array(self.facets)) )
        np.testing.assert_array_equal( sort_2D_array(np.array(polytope2.facets)), sort_2D_array(np.array(self.facets)) )

    def test_facet2vertex(self):
        polytope1 = Polytope(facets=self.facets)
        polytope2 = Polytope(facets=self.facets_with_outside)

        # useful error message?
        for polytope in (polytope1, polytope2):
            self.assertRaisesRegexp(
                                        OSError, 'No such file or directory.*nonexistentNormalizExecutable',
                                        polytope.complete_representation, normaliz='nonexistentNormalizExecutable',
                                        workdir='tmpdir_test_facet2vertex_python' + python_major_version
                                   )

        polytope1.complete_representation(workdir='tmpdir1_test_facet2vertex_python' + python_major_version)
        polytope2.complete_representation(workdir='tmpdir2_test_facet2vertex_python' + python_major_version)

        self.assertRaisesRegexp(ValueError, '(B|b)oth.*already', polytope1.complete_representation, workdir='tmpdir3_test_facet2vertex_python' + python_major_version)
        self.assertRaisesRegexp(ValueError, '(B|b)oth.*already', polytope2.complete_representation, workdir='tmpdir4_test_facet2vertex_python' + python_major_version)

        # The ordering is not important but must be fixed to compare the arrays
        np.testing.assert_array_equal( sort_2D_array(np.array(polytope1.vertices)), sort_2D_array(np.array(self.vertices)) )
        np.testing.assert_array_equal( sort_2D_array(np.array(polytope2.vertices)), sort_2D_array(np.array(self.vertices)) )
        np.testing.assert_array_equal( sort_2D_array(np.array(polytope1.facets)), sort_2D_array(np.array(self.facets)) )
        np.testing.assert_array_equal( sort_2D_array(np.array(polytope2.facets)), sort_2D_array(np.array(self.facets)) )
