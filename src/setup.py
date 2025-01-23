try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

from Cython.Build import cythonize
import numpy

setup(ext_modules=cythonize('CySPT.pyx') , compiler_directives={"language_level": "3"})
setup(ext_modules=cythonize('CyUtils.pyx'), compiler_directives={"language_level": "3"}, include_dirs=[numpy.get_include()])