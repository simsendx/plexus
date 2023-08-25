from setuptools import setup, Extension
from Cython.Build import cythonize

# To build the cython extensions, run:
# python setup.py build_ext --inplace

# Define extensions
extensions = [
    Extension("hello", ["hello.pyx"])
]

# Run setup
setup(
    ext_modules = cythonize(extensions)
)