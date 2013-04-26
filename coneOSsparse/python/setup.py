from distutils.core import setup, Extension
from glob import glob
from numpy import get_include

direct = Extension('coneOS_direct', libraries = ['m'],
                    include_dirs = ['../', get_include()],
                    define_macros = [('PYTHON', None)],
                    sources = ['coneOSmodule.c',
                        '../linAlg.c', '../cones.c', '../cs.c', 
                        '../coneOS.c', '../util.c'
                    ] + glob('../direct/*.c'),
                    extra_compile_args=['-std=c99'])

indirect = Extension('coneOS_indirect', libraries = ['m'],
                    include_dirs = ['../', get_include()],
                    define_macros = [('INDIRECT', None), ('PYTHON', None)],
                    sources = ['coneOSmodule.c',
                        '../indirect/private.c',
                        '../linAlg.c', '../cones.c', '../cs.c', 
                        '../coneOS.c', '../util.c'
                    ],
                    extra_compile_args=['-std=c99'])


setup(  name = 'coneOS',
        version = '1.0',
        description = 'This is Python package to wrap our first-order solvers',
        ext_modules = [direct, indirect],
        requires = ["numpy (>= 1.7)"])
