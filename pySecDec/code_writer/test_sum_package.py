from .sum_package import *
from ..algebra import Polynomial, ExponentiatedPolynomial
from nose.plugins.attrib import attr
import numpy as np
import sympy as sp
import unittest

#@attr('active')
class TestSumPackage(unittest.TestCase):
    #@attr('active')
    def test_coefficient_order(self):
        coefficient0 = Coefficient([],[],[1])
        np.testing.assert_array_equal(coefficient0.orders, [0])

        coefficient1 = Coefficient([Polynomial([[1,2,3],[2,1,1]],[2,1]), Polynomial([[2,2,1]],[3])],[], [2])
        np.testing.assert_array_equal(coefficient1.orders, [2])

        coefficient2 = Coefficient([],[Polynomial([[1,2,3],[2,1,1]],[2,1]), Polynomial([[2,2,1]],[3])], [2])
        np.testing.assert_array_equal(coefficient2.orders, [-2])

        coefficient3 = Coefficient([ExponentiatedPolynomial([[2,4,3],[1,2,2]],[1,1],3)],[Polynomial([[1,2,3],[2,1,2]],[2,1]), Polynomial([[2,2,1]],[3])], [2])
        np.testing.assert_array_equal(coefficient3.orders, [3])

        coefficient4 = Coefficient([ExponentiatedPolynomial([[2,4,3],[1,2,2]],[1,1],-2)],[Polynomial([[1,2,-1],[2,1,2]],[2,1]), Polynomial([[2,2,1]],[3])], [1,2])
        np.testing.assert_array_equal(coefficient4.orders, [-7,-4])
