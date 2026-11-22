from setuptools import setup
from setuptools import Extension

from Cython.Build import cythonize
import numpy

# Compiler directives belong to cythonize (passed to setup they are ignored,
# with an "Unknown distribution option" warning); include_dirs belongs to the
# Extension.
directives = {"language_level": "3"}
extensionsSPT = [Extension("CySPT", ["CySPT.pyx"])]
extensionsUtils = [Extension("CyUtils", ["CyUtils.pyx"], include_dirs=[numpy.get_include()])]

setup(ext_modules=cythonize(extensionsSPT, compiler_directives=directives))
setup(ext_modules=cythonize(extensionsUtils, compiler_directives=directives))
