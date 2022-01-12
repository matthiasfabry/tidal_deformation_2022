"""
setup script to build the cython extentions into usable python modules
do
$ python setup.py build_ext --inplace
to run the setup. requires setuptools.
"""

import glob

from setuptools import Extension, setup
import numpy

USE_CYTHON = True  # command line option, try-import, ...

ext = '.pyx' if USE_CYTHON else '.c'
files = glob.glob('c**'+ext)
extensions = list()
for file in files:
    extensions.append(Extension(file.split('.')[0], [file.split('.')[0] + ext]))

if USE_CYTHON:
    from Cython.Build import cythonize

    extensions = cythonize(extensions, annotate=True)

setup(
    ext_modules=extensions, include_dirs=[numpy.get_include()]
)
