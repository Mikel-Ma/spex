from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
import sys
import setuptools
import pybind11

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path"""
    def __str__(self):
        return pybind11.get_include()
    
class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'unix': ['-O3', '-std=c++17'],
    }

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        for ext in self.extensions:
            ext.extra_compile_args = opts
        build_ext.build_extensions(self)

ext_modules = [
    Extension(
        'spex_tequila',
        ['spex.cpp'],
        include_dirs=[
            get_pybind_include(),
        ],
        language='c++',
    ),
]

setup(
    name='spex-tequila',
    version='0.0.4',
    author='Michael Lang',
    author_email='lang-michi@t-online.de', 
    url='https://git.rz.uni-augsburg.de/qalg-a/spex', 
    description='Expectation value computation module on sparse Pauli states for Tequila, implemented in C++ using Pybind11',
    long_description=open('README.md', encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    ext_modules=ext_modules,
    cmdclass={'build_ext': build_ext},
    zip_safe=False,
    python_requires='>=3.6',
    install_requires=[
        'pybind11>=2.5.0',
        'tequila'
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: C++',
        'Operating System :: OS Independent',
        'License :: OSI Approved :: MIT License',
    ],
    packages=find_packages(),
)
