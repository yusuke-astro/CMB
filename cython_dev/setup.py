from distutils.core import setup
from Cython.Build import cythonize
import numpy as np
setup(name = 'loop_func', ext_modules = cythonize('loop_func.pyx'), include_dirs = [np.get_include()])
