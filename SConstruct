# Python packaging has layers upon layers of cruft (setuptools)
# on top of a very simple concept (zip files extracted into the
# site-packages directory). You can read [1] to get a sense of it.
#
# [1] https://blog.schuetze.link/2018/07/21/a-dive-into-packaging-native-python-extensions.html
#
# In short:
#
# 1. A source distribution (sdist) is a specially named .tar.gz
#    archive with a special 3-line file PKG-INFO; it can otherwise
#    contain anything.
#
#    This file will be extracted into a temporary location, where
#    it must be able to produce a built distribution.
#
# 2. A build distribution (bdist, aka a wheel) is a specially named
#    .whl (zip) archive with three special files:
#    * <package name>-<package version>.dist-info/WHEEL
#    * <package name>-<package version>.dist-info/METADATA
#    * <package name>-<package version>.dist-info/RECORD
#
#    This whole file will be extracted into site-packages,
#    and the dist-info/RECORD file will be used to keep
#    track of which files belong to which packages.
#
# Therefore, in this SConstruct file we must list every source
# for the sdist, and every binary needed in the bdist. Then,
# SCons will do the rest for us.

import enscons
import os
import pytoml
import subprocess

pyproject = pytoml.load(open("pyproject.toml"))

env = Environment(
    tools = ["default", "packaging", enscons.generate],
    PACKAGE_METADATA = pyproject["tool"]["enscons"],
    WHEEL_TAG = enscons.get_abi3_tag(),
    ENV = os.environ
)

def DirectoryFiles(dir):
    files = []
    for root, dirnames, filenames in os.walk(dir):
        files += File(sorted(filenames), root)
    return files

contrib = SConscript("pySecDecContrib/SConscript", exports="env")
source = DirectoryFiles("pySecDec")
sdist_extra_source = File(Split("""
    COPYING
    ChangeLog
    PKG-INFO
    README.md
    SConstruct
    pyproject.toml
"""))

if os.path.exists(".git"):
    git_id = subprocess.check_output(["git", "rev-parse", "HEAD"], encoding="utf8").strip()
    metadata_py = env.Textfile(target="pySecDec/metadata.py", source=[
        f'__version__ = version = "{pyproject["tool"]["enscons"]["version"]}"',
        f'__authors__ = authors = "{pyproject["tool"]["enscons"]["author"]}"',
        f'__commit__ = git_id = "{git_id}"'
    ])
    AlwaysBuild(metadata_py)
else:
    # We are in a giless sdist. Lets hope that metadata.py was
    # included in it.
    metadata_py = []

platformlib = env.Whl("platlib", source + contrib, root="")
bdist = env.WhlFile(source=platformlib)

# FindSourceFiles() will list every source file of every target
# defined so far.
sdist = env.SDist(source=FindSourceFiles() + metadata_py)

env.Alias("sdist", sdist)
env.Alias("bdist", bdist)
env.Alias("wheel", bdist)
env.Alias("contrib", contrib)

env.Default(contrib, sdist, bdist)
