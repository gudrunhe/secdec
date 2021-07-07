from .template_parser import *
from nose.plugins.attrib import attr
import sys, os, shutil
import unittest

python_major_version = sys.version[0]

class TestTemplateParsing(unittest.TestCase):
    def setUp(self):
        # create a temporary working directory
        i = 0
        nodir = True
        while i < 10 and nodir:
            try:
                i += 1
                self.tmpdir = 'tmpdir_test_template_parsing_%i_python'%i + python_major_version
                os.mkdir(self.tmpdir)
                nodir = False
            except OSError:
                pass
        if nodir:
            raise

    def tearDown(self):
        # delete the temporary directory
        shutil.rmtree(self.tmpdir)

    #@attr('active')
    def test_parse_template_file(self):
        # create a template file
        path_to_template_file = os.path.join(self.tmpdir, 'template')
        with open(path_to_template_file, 'w') as template_file:
            template_file.write('''
            unchanged line
            inserting integer: %(number)i
            inserting string: %(string)s''')

        # define replacements
        replacements = dict(number=1, string='Hello world')

        # parse file using `parse_template_file`
        path_to_parsed_file = os.path.join(self.tmpdir, 'parsed')
        parse_template_file(path_to_template_file, path_to_parsed_file, replacements)

        # read in parsed file
        with open(path_to_parsed_file, 'r') as parsed_file:
            parsed = parsed_file.read()

        # expected content of the parsed file
        target_parsed = '''
            unchanged line
            inserting integer: 1
            inserting string: Hello world'''

        self.assertEqual(parsed, target_parsed)

    #@attr('active')
    def test_parse_template_tree(self):
        # create template file tree
        path_to_template_tree = os.path.join(self.tmpdir, 'templates')
        os.mkdir(path_to_template_tree)
        os.mkdir(os.path.join(path_to_template_tree, 'dir1'))
        path_to_template_file_1 = os.path.join(path_to_template_tree, 'file1')
        path_to_template_file_2 = os.path.join(path_to_template_tree, 'file2')
        path_to_template_file_3 = os.path.join(path_to_template_tree, 'file3')
        path_to_template_file_4 = os.path.join(path_to_template_tree, 'dir1', 'file4')
        os.mkdir(os.path.join(path_to_template_tree, 'dir1', 'subdir1'))
        with open(path_to_template_file_1, 'w') as template_file:
            template_file.write('''
            unchanged line
            inserting integer: %(number)i
            inserting string: %(string)s''')
        with open(path_to_template_file_2, 'w') as template_file:
            template_file.write('''
            unchanged line
            inserting integer: %(number)i
            inserting string: %(string)s''')
        with open(path_to_template_file_3, 'w') as template_file:
            template_file.write('''
            unchanged line
            inserting integer: %(number)i
            inserting string: %(string)s''')
        with open(path_to_template_file_4, 'w') as template_file:
            template_file.write('''
            unchanged line
            inserting integer: %(number)i
            inserting string: %(string)s''')

        # define replacements
        replacements_in_files = dict(number=1, string='Hello world')
        replace_filenames = {
                                 'file1' : 'renamed1',
                                 'file3' : None, # do not parse
                                 'dir1' : 'renamed_dir'
                            }

        # parse file using `parse_template_file`
        path_to_parsed_tree = os.path.join(self.tmpdir, 'parsed_tree')
        path_to_parsed_file_1 = os.path.join(path_to_parsed_tree, 'renamed1')
        path_to_parsed_file_2 = os.path.join(path_to_parsed_tree, 'file2')
        path_to_parsed_file_4 = os.path.join(path_to_parsed_tree, 'renamed_dir', 'file4')
        parse_template_tree(path_to_template_tree, path_to_parsed_tree, replacements_in_files, replace_filenames)

        # should only have the parsed files (with one renamed and without 'file3') in the destination tree
        self.assertEqual(set(os.listdir(path_to_parsed_tree)), set(['renamed1', 'file2', 'renamed_dir']))
        self.assertEqual(set(os.listdir(os.path.join(path_to_parsed_tree, 'renamed_dir'))), set(['file4', 'subdir1']))

        for path_to_parsed_file in (path_to_parsed_file_1, path_to_parsed_file_2, path_to_parsed_file_4):
            print('Comparing file "' + path_to_parsed_file + '"')

            # read in parsed file
            with open(path_to_parsed_file, 'r') as parsed_file:
                parsed = parsed_file.read()

            # expected content of the first parsed file
            target_parsed = '''
            unchanged line
            inserting integer: 1
            inserting string: Hello world'''

            self.assertEqual(parsed, target_parsed)
