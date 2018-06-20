from .draw import *
from nose.plugins.attrib import attr
import sys, os
import unittest

python_major_version = sys.version[0]

class TestPlotDiagram(unittest.TestCase):
    #@attr('active')
    def test_no_powerlist(self):
        int_lines = [['m',['a',4]],['m',[4,5]],['m',['a',5]],[0,[1,2]],[0,[4,1]],[0,[2,5]]]
        ext_lines = [['p1',1],['p2',2],['p3','a']]
        try:
            plot_diagram(int_lines, ext_lines, 'tmpfile_test_plot_diagram_no_powerlist_python' + python_major_version, extension='ps')
        finally:
            os.remove('tmpfile_test_plot_diagram_no_powerlist_python' + python_major_version + '.ps')

        # input unchanged ?
        self.assertEqual(int_lines, [['m',['a',4]],['m',[4,5]],['m',['a',5]],[0,[1,2]],[0,[4,1]],[0,[2,5]]])
        self.assertEqual(ext_lines, [['p1',1],['p2',2],['p3','a']])

    #@attr('active')
    def test_with_powerlist(self):
        int_lines=[['m',[1,5]],['m',[2,6]],['M',[1,2]],['M',[3,5]],
                   ['m',[3,6]],['m',[4,6]],['M',[4,5]]]
        ext_lines=[['p1',1],['p2',2],['p3',3],['p4',4]]
        powerlist = [2,1,0,-1,-2,1,1]
        try:
            plot_diagram(int_lines, ext_lines, 'tmpfile_test_plot_diagram_with_powerlist_python' + python_major_version, powerlist, extension='ps')
        finally:
            os.remove('tmpfile_test_plot_diagram_with_powerlist_python' + python_major_version + '.ps')

        # input unchanged ?
        self.assertEqual(int_lines, [['m',[1,5]],['m',[2,6]],['M',[1,2]],['M',[3,5]],
                                     ['m',[3,6]],['m',[4,6]],['M',[4,5]]])
        self.assertEqual(ext_lines, [['p1',1],['p2',2],['p3',3],['p4',4]])
