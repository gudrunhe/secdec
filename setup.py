from __future__ import print_function

# bootstrap: download setuptools 3.3 if needed
from ez_setup import use_setuptools
use_setuptools()

from setuptools import setup, find_packages

# get the version number and author list from package
import sys
sys.path.insert(0, 'pySecDec')
from metadata import __version__, __authors__

setup(
    name='pySecDec',
    packages=find_packages(),
    package_data={'pySecDec.code_writer': ['templates' + '/*' * depth for depth in range(1,10)]},
    version=__version__,
    author=__authors__,
    author_email='secdec@projects.hepforge.org',
    url='secdec.hepforge.org',
    description='An implementation of "Sector Decomposition" (see arXiv:0803.4177).',
    # long_description=long_description, #TODO: write an abstract of this package
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
