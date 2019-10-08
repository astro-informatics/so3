import sys
import os
import shutil

from distutils.core import setup, Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

import numpy

# clean previous build
for root, dirs, files in os.walk('./src/python/', topdown=False):
    for name in dirs:
        if (name == 'build'):
            shutil.rmtree(name)

include_dirs = [
    numpy.get_include(),
    './include/c',
    os.environ['SSHT'] + '/src/c',
]

extra_link_args = [
    '-L./build',
    '-L' + os.environ['FFTW'] + '/lib',
    # depending whether used cmake
    '-L' + os.environ['SSHT'] + '/build',
    # or make
    '-L' + os.environ['SSHT'] + '/lib/c',
]

setup(
    name='pyso3',
    version='2.0',
    cmdclass={'build_ext': build_ext},
    ext_modules=cythonize([Extension(
        'src/python/pyso3',
        sources=['src/python/pyso3.pyx'],
        include_dirs=include_dirs,
        libraries=['so3', 'ssht', 'fftw3'],
        extra_link_args=extra_link_args,
        extra_compile_args=[]
    )])
)