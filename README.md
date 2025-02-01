# SPEX
spex is an expectation value computation module on sparse Pauli states for Tequila, implemented in C++ using Pybind11. It provides computation of expectation values, inner products, and application of exponential Pauli operators on sparse quantum states.

# Install on MacOS
the default c++ compiler (clang) does not support OpenMP. You can however install spex by using GNU compilers, one way of getting them:

```bash
brew gcc@14
```

Change into the directory where `setup.py` is located and install with 

```bash
CC=gcc-14 CXX=g++-14 pip install .
```

# Install on Linux

Assuming that c++ compilers are on the system
```bash
pip install .
```