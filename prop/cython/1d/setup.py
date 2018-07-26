from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np 

setup(
    
    ext_modules=[
        Extension('prop1d_cython',
        sources=['prop1d_cython.pyx'],
        extra_compile_args=['-O3','-fopenmp','-mavx'],
        language='c')
        ],
    
    include_dirs = [np.get_include()],
    
    cmdclass = {'build_ext': build_ext}
    
)
