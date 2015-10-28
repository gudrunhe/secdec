"""Unit tests for the miscellaneous routines"""

from .misc import *
import unittest

class TestMisc(unittest.TestCase):
    def test_powerset_range(self):
        self.assertEqual(list(powerset(range(0))), [()])
        self.assertEqual(list(powerset(range(1))), [(),(0,)])
        self.assertEqual(list(powerset(range(2))), [(),(0,),(1,),(0,1)])
        self.assertEqual(list(powerset(range(3))), [(),(0,),(1,),(2,),(0,1),(0,2),(1,2),(0,1,2)])
        self.assertEqual(list(powerset(range(4))), [(),(0,),(1,),(2,),(3,),(0,1),(0,2),(0,3),(1,2),(1,3),(2,3),(0,1,2),(0,1,3),(0,2,3),(1,2,3),(0,1,2,3)])
