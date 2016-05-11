# Use the following command to compile:
# python cython_psd_make.py build_ext --inplace
from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("cython_psd.pyx")
)