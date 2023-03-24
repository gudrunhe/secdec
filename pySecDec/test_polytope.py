from .polytope import *
from .misc import argsort_2D_array
import sys
import unittest
import pytest

python_major_version = sys.version[0]

def sort_2D_array(array):
    'Use the .misc.argsort_2D_array function to sort an array'
    return array[argsort_2D_array(array)]

#@pytest.mark.active
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
        with pytest.raises(TypeError, match='either.*vertices.*or.*facets'):
            Polytope(facets=self.facets, vertices=self.vertices)
        with pytest.raises(TypeError, match='either.*vertices.*or.*facets'):
            Polytope()

        inarray = np.array([[1,2,3],
                            [4,5,6]])

        polytope1 = Polytope(facets=inarray)
        polytope2 = Polytope(vertices=inarray)

        inarray[0,0] = 10

        assert polytope1.facets[0,0] == 1
        assert polytope2.vertices[0,0] == 1

    #@pytest.mark.active
    @pytest.mark.slow
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

    #@pytest.mark.active
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
            with pytest.raises(AssertionError, match='complete_representation.*first'):
                polytope.vertex_incidence_lists()

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
            assert len(incidence.keys()) == 6
            for key, value in incidence.items():
                assert key in target_incidence_lists.keys()

                # The ordering is not important but must be fixed to compare the arrays
                np.testing.assert_array_equal(sort_2D_array(polytope.facets[value]), sort_2D_array(target_incidence_lists[key]))

    #@pytest.mark.active
    def test_vertex2facet(self):
        polytope1 = Polytope(vertices=self.vertices)
        polytope2 = Polytope(vertices=self.vertices_with_inside)

        # useful error message?
        for polytope in (polytope1, polytope2):
            with pytest.raises(OSError, match='No such file or directory.*nonexistentNormalizExecutable'):
                polytope.complete_representation(normaliz='nonexistentNormalizExecutable',
                                        workdir='tmpdir_test_vertex2facet_python' + python_major_version)

        polytope1.complete_representation(workdir='tmpdir1_test_vertex2facet_python' + python_major_version)
        polytope2.complete_representation(workdir='tmpdir2_test_vertex2facet_python' + python_major_version)

        with pytest.raises(ValueError, match='(B|b)oth.*already'):
            polytope1.complete_representation(workdir='tmpdir3_test_vertex2facet_python' + python_major_version)
        with pytest.raises(ValueError, match='(B|b)oth.*already'):
            polytope2.complete_representation(workdir='tmpdir4_test_vertex2facet_python' + python_major_version)

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
            with pytest.raises(OSError, match='No such file or directory.*nonexistentNormalizExecutable'):
                polytope.complete_representation(normaliz='nonexistentNormalizExecutable',
                                        workdir='tmpdir_test_facet2vertex_python' + python_major_version)

        polytope1.complete_representation(workdir='tmpdir1_test_facet2vertex_python' + python_major_version)
        polytope2.complete_representation(workdir='tmpdir2_test_facet2vertex_python' + python_major_version)

        with pytest.raises(ValueError, match='(B|b)oth.*already'):
            polytope1.complete_representation(workdir='tmpdir3_test_facet2vertex_python' + python_major_version)
        with pytest.raises(ValueError, match='(B|b)oth.*already'):
            polytope2.complete_representation(workdir='tmpdir4_test_facet2vertex_python' + python_major_version)

        # The ordering is not important but must be fixed to compare the arrays
        np.testing.assert_array_equal( sort_2D_array(np.array(polytope1.vertices)), sort_2D_array(np.array(self.vertices)) )
        np.testing.assert_array_equal( sort_2D_array(np.array(polytope2.vertices)), sort_2D_array(np.array(self.vertices)) )
        np.testing.assert_array_equal( sort_2D_array(np.array(polytope1.facets)), sort_2D_array(np.array(self.facets)) )
        np.testing.assert_array_equal( sort_2D_array(np.array(polytope2.facets)), sort_2D_array(np.array(self.facets)) )

    #@pytest.mark.active
    def test_equations(self):
        polytope1 = Polytope(vertices=[[0,1],[1,0]])
        polytope2 = Polytope(vertices=[[0,1],[1,0],[1,1]])
        polytope3 = Polytope(facets=[[0,1,0],[0,-1,0],[1,0,0],[-1,0,1]])

        polytope1.complete_representation(workdir='tmpdir_test_equations1_python' + python_major_version)
        polytope2.complete_representation(workdir='tmpdir_test_equations2_python' + python_major_version)
        polytope3.complete_representation(workdir='tmpdir_test_equations3_python' + python_major_version)

        np.testing.assert_array_equal(polytope1.equations, np.array([[1,1,-1]]))
        np.testing.assert_array_equal(polytope2.equations, np.array([]))
        np.testing.assert_array_equal(polytope3.equations, np.array([[0,1,0]]))
        np.testing.assert_array_equal(sort_2D_array(polytope3.facets), sort_2D_array(np.array([[-1,0,1],[1,0,0]])))


    #@pytest.mark.active
    def test_triangulate(self):
        # basic consistency checks working?
        simplicial_cone = [[ 1,  0,  0], [ 0,  1,  0], [ 0, -1, -1]]
        with pytest.raises(ValueError, match='simplicial.*already'):
            triangulate(simplicial_cone)
        wrong_dimensionality = [ 1,  0,  0]
        with pytest.raises(AssertionError, match='(M|m)ust.*two.*dim'):
            triangulate(wrong_dimensionality)
        two_rays = [[ 1,  0,  0], [ 0,  1,  0]]
        with pytest.raises(AssertionError, match='(M|m)ust.*at least.*dim'):
            triangulate(two_rays)


        cone = [[ 1,  0,  0], [ 0,  1,  0], [ 0, -1, -1], [-1,  0, -1]]
        cone_normal = [[ -1, 1, 1], [ 1, 0, 0], [ 0, 1, 0], [ 0, 0, 1]]

        # useful error message?
        with pytest.raises(OSError, match='No such file or directory.*nonexistentNormalizExecutable'):
            triangulate(cone, normaliz='nonexistentNormalizExecutable',
                                    workdir='tmpdir_test_triangulate_python' + python_major_version)

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

    def test_triangulate_0D(self):
        cone = np.array([[1,0],[-1,0],[0,1],[0,-1]])
        target_cone = np.array([])
        triangulated_cone = triangulate(cone, workdir='tmpdir_test_triangulate_0D_python' + python_major_version, switch_representation=True)
        np.testing.assert_array_equal(triangulated_cone, target_cone)
