from __future__ import print_function

# bootstrap: download setuptools 3.3 if needed
from ez_setup import use_setuptools
use_setuptools()

from setuptools import setup, find_packages

# get the version number and author list from package
import sys
sys.path.insert(0, 'pySecDec')
from metadata import __version__, __authors__

# collect the c++ template files needed by `pySecDec.code_writer.make_package`
import os
cpp_template_files = []
relpath_to_code_writer_module = os.path.join('pySecDec','code_writer')
for dir_name, subdir_list, file_list in os.walk(relpath_to_code_writer_module):
    relpath_to_templates = os.path.relpath(dir_name, relpath_to_code_writer_module)
    for data_file in file_list:
        cpp_template_files.append( os.path.join(relpath_to_templates, data_file) )

setup(
    name='pySecDec',
    packages=find_packages(),
    package_data={'pySecDec.code_writer': cpp_template_files},
    version=__version__,
    author=__authors__,
    author_email='secdec@projects.hepforge.org',
    url='secdec.hepforge.org',
    description='An implementation of "Sector Decomposition" (see arXiv:hep-ph/0004013, arXiv:0803.4177).', # TODO: cite paper for this SecDec version
    license='GPLv3',
    install_requires=['numpy>=1.6', 'sympy', 'setuptools>=3.3'],
    extras_require={'testing': ['nose'], 'documentation': ['sphinx']},
    classifiers=['Development Status :: 3 - Alpha',
                 'Environment :: Console',
                 'Intended Audience :: Science/Research',
                 'Natural Language :: English',
                 'Operating System :: Unix',
                 'Programming Language :: Python :: 2.7',
                 'Programming Language :: Python :: 3',
                 'Programming Language :: C++',
                 'Topic :: Scientific/Engineering :: Physics',
                 ],
    platforms=['Unix'],
    )
