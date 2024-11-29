from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
    Pybind11Extension(
        'spex_tequila',
        ['spex.cpp'],
    ),
]

setup(
    name='spex-tequila',
    version='0.1',
    author='Michael Lang',
    description='Expectation value computation module on sparse pauli states for tequila, implemented in C++ using Pybind11',
    ext_modules=ext_modules,
    cmdclass={'build_ext': build_ext},
    zip_safe=False,
)
