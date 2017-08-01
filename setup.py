from __future__ import print_function

# bootstrap: download setuptools 3.3 if needed
from ez_setup import use_setuptools
use_setuptools()

from setuptools import setup, find_packages

# collect the c++ template files needed by `pySecDec.code_writer.make_package`
import os
cpp_template_files = []
relpath_to_code_writer_module = os.path.join('pySecDec','code_writer')
for dir_name, subdir_list, file_list in os.walk(relpath_to_code_writer_module):
    relpath_to_templates = os.path.relpath(dir_name, relpath_to_code_writer_module)
    for data_file in file_list:
        cpp_template_files.append( os.path.join(relpath_to_templates, data_file) )

# prevent python's distutils to try (and potentially fail) hard linking
# this is due to a bug in python's distutils, see: http://bugs.python.org/issue8876
del os.link

# get the git commit id
from subprocess import check_output as shell_call, CalledProcessError, STDOUT
try:
    git_id = shell_call(['git','rev-parse','HEAD'], cwd=os.path.split(os.path.abspath(__file__))[0], stderr=STDOUT).strip()
except CalledProcessError:
    git_id = 'UNKNOWN'
# in python3, the return type of `shell_call` may be `bytes` but we need `str`
if not isinstance(git_id, str):
    git_id = git_id.decode()

# hard-code the current git commit id into `metadata.py`
import shutil
shutil.copy(os.path.join('pySecDec','metadata.py'), os.path.join('pySecDec','metadata.py.orig'))
try:
    # hard-code the current git commit id into `metadata.py` (continued)
    with open(os.path.join('pySecDec','metadata.py'), 'r') as metadata_file:
        metadata = metadata_file.read()
    metadata = metadata.replace(
"""
from subprocess import check_output as shell_call
import os
# The following line is replaced by `setup.py` --> hard-code the commit id in distributions
git_id = shell_call(['git','rev-parse','HEAD'], cwd=os.path.split(os.path.abspath(__file__))[0]).strip()
# in python3, the return type of `shell_call` may be `bytes` but we need `str`
if not isinstance(git_id, str):
    git_id = git_id.decode()""", "git_id = '" + git_id + "'")
    with open(os.path.join('pySecDec','metadata.py'), 'w') as metadata_file:
        metadata_file.write(metadata)

    # get the version number and author list from package
    import sys
    sys.path.insert(0, 'pySecDec')
    from metadata import __version__, __authors__

    setup(
        name='pySecDec',
        packages=find_packages(),
        package_data={'pySecDec.code_writer': cpp_template_files},
        version=__version__,
        author=__authors__,
        author_email='secdec@projects.hepforge.org',
        url='secdec.hepforge.org',
        description='An implementation of "Sector Decomposition" (see arXiv:1703.09692, arXiv:hep-ph/0004013, arXiv:0803.4177).',
        license='GPLv3',
        install_requires=['numpy>=1.6', 'sympy>=0.7.6', 'setuptools>=3.3'],
        extras_require={'testing': ['nose'], 'documentation': ['sphinx>=1.6.3']},
        classifiers=['Development Status :: 5 - Production/Stable',
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
finally:
    # restore the original `metadata.py`
    shutil.move(os.path.join('pySecDec','metadata.py.orig'), os.path.join('pySecDec','metadata.py'))
