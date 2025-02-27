try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

from Cython.Build import cythonize
import numpy

extensionsSPT = [Extension("CySPT", ["CySPT.pyx"])]
extensionsUtils = [Extension("CyUtils", ["CyUtils.pyx"])]

setup(ext_modules=cythonize(extensionsSPT) , compiler_directives={"language_level": "3"})
setup(ext_modules=cythonize(extensionsUtils), compiler_directives={"language_level": "3"}, include_dirs=[numpy.get_include()])
