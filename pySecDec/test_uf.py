"""Unit tests for the U, F routines"""

from .uf import *
import sympy as sp
import unittest

class TestUF(unittest.TestCase):
    def test_box_1l(self):
        box_1l = uf(['k1'],['k1**2','(k1-p1)**2'])
        zerou = sp.sympify("""-x0 - x1""") - box_1l[0]
        zerof = sp.sympify("""-p1**2*x0*x1""") - box_1l[1]
        zerou.simplify()
        zerof.simplify()
        self.assertEqual(zerou,0)
        self.assertEqual(zerof,0)

    def test_box_2l(self):
        box_2l = uf(['k1','k2'],
                    ['k1**2','(k1+p2)**2','(k1-p1)**2','(k1-k2)**2',
                    '(k2+p2)**2','(k2-p1)**2','(k2+p2+p3)**2'])
        zerou = sp.sympify("""x0*x3 + x0*x4 + x0*x5 + x0*x6 + x1*x3 +
                              x1*x4 + x1*x5 + x1*x6 + x2*x3 + x2*x4 +
                              x2*x5 + x2*x6 + x3*x4 + x3*x5 + x3*x6""") - box_2l[0]
        zerof = sp.sympify("""p1**2*x0*x2*x3 + p1**2*x0*x2*x4 + p1**2*x0*x2*x5 +
                              p1**2*x0*x2*x6 + p1**2*x0*x3*x5 + p1**2*x0*x4*x5 +
                              p1**2*x0*x5*x6 + p1**2*x1*x2*x3 + p1**2*x1*x2*x4 +
                              p1**2*x1*x2*x5 + p1**2*x1*x2*x6 + p1**2*x1*x3*x5 +
                              p1**2*x1*x4*x5 + p1**2*x1*x5*x6 + p1**2*x2*x3*x4 +
                              p1**2*x2*x3*x6 + p1**2*x2*x4*x5 + p1**2*x2*x5*x6 +
                              p1**2*x3*x4*x5 + p1**2*x3*x5*x6 + 2*p1*p2*x0*x4*x5 +
                              2*p1*p2*x0*x5*x6 + 2*p1*p2*x1*x2*x3 +
                              2*p1*p2*x1*x2*x4 + 2*p1*p2*x1*x2*x5 +
                              2*p1*p2*x1*x2*x6 + 2*p1*p2*x1*x3*x5 +
                              2*p1*p2*x1*x4*x5 + 2*p1*p2*x1*x5*x6 +
                              2*p1*p2*x2*x3*x4 + 2*p1*p2*x2*x3*x6 +
                              2*p1*p2*x2*x4*x5 + 2*p1*p2*x2*x5*x6 +
                              2*p1*p2*x3*x4*x5 + 2*p1*p2*x3*x5*x6 +
                              2*p1*p3*x0*x5*x6 + 2*p1*p3*x1*x5*x6 +
                              2*p1*p3*x2*x3*x6 + 2*p1*p3*x2*x5*x6 +
                              2*p1*p3*x3*x5*x6 + p2**2*x0*x1*x3 + p2**2*x0*x1*x4 +
                              p2**2*x0*x1*x5 + p2**2*x0*x1*x6 + p2**2*x0*x3*x4 +
                              p2**2*x0*x3*x6 + p2**2*x0*x4*x5 + p2**2*x0*x5*x6 +
                              p2**2*x1*x2*x3 + p2**2*x1*x2*x4 + p2**2*x1*x2*x5 +
                              p2**2*x1*x2*x6 + p2**2*x1*x3*x5 + p2**2*x1*x4*x5 +
                              p2**2*x1*x5*x6 + p2**2*x2*x3*x4 + p2**2*x2*x3*x6 +
                              p2**2*x2*x4*x5 + p2**2*x2*x5*x6 + p2**2*x3*x4*x5 +
                              p2**2*x3*x5*x6 + 2*p2*p3*x0*x3*x6 +
                              2*p2*p3*x0*x5*x6 + 2*p2*p3*x1*x5*x6 +
                              2*p2*p3*x2*x3*x6 + 2*p2*p3*x2*x5*x6 +
                              2*p2*p3*x3*x5*x6 + p3**2*x0*x3*x6 + p3**2*x0*x4*x6 +
                              p3**2*x0*x5*x6 + p3**2*x1*x3*x6 + p3**2*x1*x4*x6 +
                              p3**2*x1*x5*x6 + p3**2*x2*x3*x6 + p3**2*x2*x4*x6 +
                              p3**2*x2*x5*x6 + p3**2*x3*x4*x6 + p3**2*x3*x5*x6""") - box_2l[1]
        zerou.simplify()
        zerof.simplify()
        self.assertEqual(zerou,0)
        self.assertEqual(zerof,0)
