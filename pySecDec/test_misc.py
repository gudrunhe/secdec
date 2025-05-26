from .misc import *
from .algebra import Polynomial
from multiprocessing import Pool
import numpy as np
import unittest
import pytest

#@pytest.mark.active
class TestMakeCPPList(unittest.TestCase):
    #@pytest.mark.active
    def test_make_cpp_list(self):
        python_list = ['a', 'b', 'c']
        cpp_list = make_cpp_list(python_list)
        target_cpp_list = '"a","b","c"'
        assert cpp_list == target_cpp_list

    #@pytest.mark.active
    def test_make_cpp_list_empty(self):
        python_list = []
        cpp_list = make_cpp_list(python_list)
        target_cpp_list = ""  # empty string
        assert cpp_list == target_cpp_list


class TestSort(unittest.TestCase):
    def test_sort_2D_array(self):
        in_array = [[1,2,3],
                    [2,3,4],
                    [1,2,3]]
        target_sorted_array = [[1,2,3],
                               [1,2,3],
                               [2,3,4]]
        calculated_sort_indices = argsort_2D_array(in_array)

        in_array = np.asarray(in_array)
        np.testing.assert_array_equal(in_array[calculated_sort_indices], target_sorted_array)

    #@pytest.mark.active
    def test_sort_ND(self):
        in_array =          [[[0, 2],
                              [1, 0]],

                             [[1, 1],
                              [0, 1]],

                             [[0, 1],
                              [1, 1]],

                             [[1, 1],
                              [0, 1]],

                             [[0, 1],
                              [1, 1]]]

        target_sorted_array = np.array([[[0, 1],
                                         [1, 1]],

                                        [[0, 1],
                                         [1, 1]],

                                        [[0, 2],
                                         [1, 0]],

                                        [[1, 1],
                                         [0, 1]],

                                        [[1, 1],
                                         [0, 1]]])

        calculated_sort_indices = argsort_ND_array(in_array)
        np.testing.assert_array_equal(np.array(in_array)[calculated_sort_indices], target_sorted_array)

    #@pytest.mark.active
    def test_error_message(self):
        with pytest.raises(AssertionError, match='array.*two dimensional'):
            argsort_2D_array([1,2,3,4])
        with pytest.raises(AssertionError, match='array.*two or higher dimensional'):
            argsort_ND_array([1,2,3,4])

class TestPowerset(unittest.TestCase):
    def test_powerset_range(self):
        assert list(powerset(range(0))) == [()]
        assert list(powerset(range(1))) == [(),(0,)]
        assert list(powerset(range(2))) == [(),(0,),(1,),(0,1)]
        assert list(powerset(range(3))) == [(),(0,),(1,),(2,),(0,1),(0,2),(1,2),(0,1,2)]
        assert list(powerset(range(4))) == [(),(0,),(1,),(2,),(3,),(0,1),(0,2),(0,3),(1,2),(1,3),(2,3),(0,1,2),(0,1,3),(0,2,3),(1,2,3),(0,1,2,3)]

    def test_powerset_without_empty(self):
        assert list(powerset(range(0),min_length=1)) == []
        assert list(powerset(range(1),min_length=1)) == [(0,)]
        assert list(powerset(range(2),min_length=1)) == [(0,),(1,),(0,1)]
        assert list(powerset(range(3),min_length=1)) == [(0,),(1,),(2,),(0,1),(0,2),(1,2),(0,1,2)]
        assert list(powerset(range(4),min_length=1)) == [(0,),(1,),(2,),(3,),(0,1),(0,2),(0,3),(1,2),(1,3),(2,3),(0,1,2),(0,1,3),(0,2,3),(1,2,3),(0,1,2,3)]

    #@pytest.mark.active
    def test_strided_powerset(self):
        assert list(powerset(range(4),stride=2)) == [(),(0,1),(0,2),(0,3),(1,2),(1,3),(2,3),(0,1,2,3)]
        assert list(powerset(range(4),stride=3)) == [(),(0,1,2),(0,1,3),(0,2,3),(1,2,3)]
        assert list(powerset(range(4),stride=2,min_length=2)) == [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3),(0,1,2,3)]

    def test_minimal_powerset(self):
        assert list(powerset(range(4),min_length=1)) == [(0,),(1,),(2,),(3,),(0,1),(0,2),(0,3),(1,2),(1,3),(2,3),(0,1,2),(0,1,3),(0,2,3),(1,2,3),(0,1,2,3)]
        assert list(powerset(range(4),min_length=2)) == [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3),(0,1,2),(0,1,3),(0,2,3),(1,2,3),(0,1,2,3)]
        assert list(powerset(range(4),stride=2,min_length=3)) == [(0,1,2),(0,1,3),(0,2,3),(1,2,3)]

class TestRangecomb(unittest.TestCase):
    #@pytest.mark.active
    def test_error_messages(self):
        with pytest.raises(TypeError):
            rangecomb(['a','b','c'], [1,2,3])
        with pytest.raises(TypeError):
            rangecomb([1,2,3], ['a','b','c'])
        with pytest.raises(AssertionError, match='low.*vector'):
            rangecomb([[-1,-2],[1,2]], [-1,-2])
        with pytest.raises(AssertionError, match='high.*vector'):
            rangecomb([-1,-2], [[-1,-2],[1,2]])
        with pytest.raises(AssertionError, match='low.*high.*length'):
            rangecomb([-1,-2], [0])

    #@pytest.mark.active
    def test_3_dimensions(self):
        low = [-2,-2, 1]
        high = [0, 1, 2]

        orders = list( rangecomb(low, high) )
        target_orders = [(-2, -2, 1), (-2, -2, 2), (-2, -1, 1),
                         (-2, -1, 2), (-2, 0, 1), (-2, 0, 2),
                         (-2, 1, 1), (-2, 1, 2), (-1, -2, 1),
                         (-1, -2, 2), (-1, -1, 1), (-1, -1, 2),
                         (-1, 0, 1), (-1, 0, 2), (-1, 1, 1),
                         (-1, 1, 2), (0, -2, 1), (0, -2, 2),
                         (0, -1, 1), (0, -1, 2), (0, 0, 1),
                         (0, 0, 2), (0, 1, 1), (0, 1, 2)]

        assert orders == target_orders

class TestMissing(unittest.TestCase):
    #@pytest.mark.active
    def test_missing(self):
        assert missing([1,2,3], [1]) == [2,3]
        assert missing([1,2,3], [1,2]) == [3]
        assert missing([1,2,3,1], [1,2]) == [3,1]
        assert missing(['a','b','c','d','e','f'], ['a','e','d']) == ['b','c','f']
        with pytest.raises(ValueError):
            missing([1,2,3], [1,'a'])

    #@pytest.mark.active
    def test_in_combination_with_powerset(self):
        full_set = list(range(4))
        powerset_4_generator = powerset(full_set)
        target_powerset_4 = [(),(0,),(1,),(2,),(3,),(0,1),(0,2),(0,3),(1,2),(1,3),(2,3),(0,1,2),(0,1,3),(0,2,3),(1,2,3),(0,1,2,3)]

        target_missing_in_powerset_4 = list(target_powerset_4)
        target_missing_in_powerset_4.reverse()

        for i, powerset_item in enumerate(powerset_4_generator):
            print(i)
            assert powerset_item == target_powerset_4[i]
            assert missing(full_set, powerset_item) == list(target_missing_in_powerset_4[i])

class TestAllPairs(unittest.TestCase):
    #@pytest.mark.active
    def test_error_input_not_even(self):
        with pytest.raises(AssertionError, match='even'):
            all_pairs([1,2,3])

    #@pytest.mark.active
    def test_partitioning(self):
        # list input
        lst = [1,2,3,4]
        assert list(all_pairs(lst)) == [[(1,2),(3,4)], [(1,3),(2,4)], [(1,4),(2,3)]]

        # iterator input
        generator = (i for i in lst)
        assert list(all_pairs(generator)) == [[(1,2),(3,4)], [(1,3),(2,4)], [(1,4),(2,3)]]

#@pytest.mark.active
class TestDet(unittest.TestCase):
    def check_calculation(self, parallel):
        if parallel:
            f = lambda M: parallel_det(M, Pool())
        else:
            f = det

        M1 = [[1,1],
              [2,2]]
        target_det_M1 = 0
        assert f(M1) == target_det_M1

        M2 = [[1,2,3],
              [0,2,2],
              [0,0,2]]
        target_det_M2 = 4
        assert f(M2) == target_det_M2

    def test_sequential(self):
        self.check_calculation(parallel=False)

    def test_parallel(self):
        self.check_calculation(parallel=True)

    def test_error_messages(self):
        # wrong shape
        M1 = [1,2,3]
        M2 = [[[1,2],[3,4]],[[5,6],[7,8]],[[9,10],[11,12]]]
        for M in (M1,M2):
            with pytest.raises(AssertionError, match='M.*two dimensional'):
                det(M)

        # not square
        M3 = [[1,2,3],
              [0,2,2]]
        with pytest.raises(AssertionError, match='M.*square'):
            det(M3)

    def test_1x1_special_case(self):
        M = [[1]]
        assert det(M) == 1
        assert np.issubdtype(type(det(M)), np.array(M).dtype)

class TestAdjugate(unittest.TestCase):
    def test_calculation(self):
        M = [[1,2,3],
             [4,5,6],
             [7,8,9]]
        target_adj_M = [[-3,   6, -3],
                        [ 6, -12,  6],
                        [-3,   6, -3]]
        np.testing.assert_array_equal(adjugate(M), target_adj_M)

    def test_error_messages(self):
        # wrong shape
        M1 = [1,2,3]
        M2 = [[[1,2],[3,4]],[[5,6],[7,8]],[[9,10],[11,12]]]
        for M in (M1,M2):
            with pytest.raises(AssertionError, match='M.*two dimensional'):
                adjugate(M)

        # not square
        M3 = [[1,2,3],
              [0,2,2]]
        with pytest.raises(AssertionError, match='M.*square'):
            adjugate(M3)

    def test_1x1_special_case(self):
        M1 = [[1 ]]
        M2 = [[5 ]]
        M3 = [[2.]]
        for M in (M1,M2,M3):
            adj_M = adjugate(M)
            assert adj_M == 1
            assert isinstance(adj_M, np.ndarray)

        # should keep data type of input
        assert np.issubdtype(adjugate(M1).dtype, np.array(M1).dtype)
        assert np.issubdtype(adjugate(M2).dtype, np.array(M2).dtype)
        assert np.issubdtype(adjugate(M3).dtype, np.array(M3).dtype)

class TestCachedProperty(unittest.TestCase):
    #@pytest.mark.active
    def test_usability(self):
        class C(object):
            'Sum up the numbers from one to `N`.'
            def __init__(self, N):
                self.N = N
            @cached_property
            def sum(self):
                result = 0
                for i in range(1, self.N + 1):
                    result += i
                return result
        to_ten = C(10)
        assert to_ten.sum == 1+2+3+4+5+6+7+8+9+10

    #@pytest.mark.active
    def test_called_only_once(self):
        class C(object):
            def __init__(self):
                self.counter = 0
            @cached_property
            def prop(self):
                self.counter += 1
                return 5 + 5
        instance = C()

        assert instance.counter == 0
        for i in range(10):
            # the method should be called only once
            assert instance.prop == 10
            assert instance.counter == 1

class TestFlatten(unittest.TestCase):
    def setUp(self):
        self.p0 = Polynomial([[-1,0], [0,0]], ['A', 'B'])
        self.p1 = Polynomial([[0,-2], [0,-1], [0,0]], [self.p0, self.p0, self.p0])
        self.p2 = Polynomial([[0,-2], [0,-1], [0,0]], [self.p0, self.p0, 'B'])

        self.target_flattened_polynomial_1 = Polynomial([[-1,-2], [0,-2],
                                                         [-1,-1], [0,-1],
                                                         [-1, 0], [0, 0]],

                                                         ['A','B']*3)


        self.target_flattened_polynomial_2 = Polynomial([[-1,-2], [0,-2],
                                                         [-1,-1], [0,-1],
                                                                  [0, 0]],

                                                         ['A','B']*2+['B'])


        self.pmix_0 = Polynomial([[ 1,2], [ 0,0]], ['A', 'B'])
        self.pmix   = Polynomial([[-1,0], [-1,0]], [self.pmix_0, 'B'])

        self.target_flattened_pmix =  Polynomial([[0,2], [-1,0], [-1,0]], ['A','B','B'])

    #@pytest.mark.active
    def test_plain_with_depth(self):
        flattened_p1 = flatten(self.p1, 1)
        np.testing.assert_array_equal(flattened_p1.expolist, self.target_flattened_polynomial_1.expolist)
        np.testing.assert_array_equal(flattened_p1.coeffs, self.target_flattened_polynomial_1.coeffs)

    #@pytest.mark.active
    def test_plain_without_depth(self):
        flattened_p1 = flatten(self.p1)
        np.testing.assert_array_equal(flattened_p1.expolist, self.target_flattened_polynomial_1.expolist)
        np.testing.assert_array_equal(flattened_p1.coeffs, self.target_flattened_polynomial_1.coeffs)

    #@pytest.mark.active
    def test_mixed_with_depth(self):
        flattened_p2 = flatten(self.p2, 1)
        np.testing.assert_array_equal(flattened_p2.expolist, self.target_flattened_polynomial_2.expolist)
        np.testing.assert_array_equal(flattened_p2.coeffs, self.target_flattened_polynomial_2.coeffs)

    #@pytest.mark.active
    def test_mixed_without_depth(self):
        flattened_p2 = flatten(self.p2)
        np.testing.assert_array_equal(flattened_p2.expolist, self.target_flattened_polynomial_2.expolist)
        np.testing.assert_array_equal(flattened_p2.coeffs, self.target_flattened_polynomial_2.coeffs)

    #@pytest.mark.active
    def test_multi_mixed_with_depth(self):
        flattened_pmix = flatten(self.pmix, 1)
        np.testing.assert_array_equal(flattened_pmix.expolist, self.target_flattened_pmix.expolist)
        np.testing.assert_array_equal(flattened_pmix.coeffs, self.target_flattened_pmix.coeffs)

    #@pytest.mark.active
    def test_multi_mixed_without_depth(self):
        flattened_pmix = flatten(self.pmix)
        np.testing.assert_array_equal(flattened_pmix.expolist, self.target_flattened_pmix.expolist)
        np.testing.assert_array_equal(flattened_pmix.coeffs, self.target_flattened_pmix.coeffs)

class TestPrefactorExpansion(unittest.TestCase):
    #@pytest.mark.active
    def test_error_message(self):
        with pytest.raises(AssertionError, match='variable.*symbol'):
            lowest_order('a+b', 'a+b')

    #@pytest.mark.active
    def test_lowest_order_one_variable_constant(self):
        expression = '''
                            -1 * 2**6 * exp(EulerGamma)**(-2*eps) / sqrt(pi) * gamma(-4*eps)
                            / gamma(1-eps) / gamma(1/2-eps) * 4 * 2**(-1-2*eps) * 2 * 16
                            / pi / 2 * (-1) * eps * 2**(-4*eps) * 4**(1+eps) * 64
                     '''
        expression_lowest_eps_order = lowest_order(expression, 'eps')
        assert expression_lowest_eps_order == 0 # no pole, expansion starts at constant order

    #@pytest.mark.active
    def test_lowest_order_one_variable_simple_pole(self):
        expression = '''
                            gamma(-5*eps - 2) * exp(EulerGamma * eps)
                     '''
        expression_lowest_eps_order = lowest_order(expression, 'eps')
        assert expression_lowest_eps_order == -1 # 1/eps pole from gamma function

    #@pytest.mark.active
    def test_lowest_order_one_variable_double_pole(self):
        expression = '''
                            gamma(-5*eps - 2)**2 * exp(EulerGamma * eps)
                     '''
        expression_lowest_eps_order = lowest_order(expression, 'eps')
        assert expression_lowest_eps_order == -2 # 1/eps**2 pole from the squared gamma function

    #@pytest.mark.active
    def test_lowest_order_one_variable_no_constant_order(self):
        expression = '''
                            gamma(-5*eps - 2) * exp(EulerGamma * eps) * (eps**3 + 15 * eps**2)
                     '''
        expression_lowest_eps_order = lowest_order(expression, 'eps')
        assert expression_lowest_eps_order == +1 # starts at ``eps**1`` order
