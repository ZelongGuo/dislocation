<!-- ![PyPI - Python Version](https://img.shields.io/pypi/pyversions/dislocs) -->
<!-- ![PyPI Version](https://img.shields.io/pypi/v/dislocs) -->
<!-- ![pypi downloads](https://img.shields.io/pypi/dm/dislocs) -->
<!-- ![License](https://img.shields.io/github/license/ZelongGuo/dislocs) -->
<!-- ![Last Commit](https://img.shields.io/github/last-commit/ZelongGuo/dislocs) -->

<!-- [![Python Versions](https://img.shields.io/pypi/pyversions/dislocs?logo=python&logoColor=white)](https://pypi.org/project/dislocs/) -->
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/dislocs?logo=python&logoColor=white)
[![PyPI Version](https://img.shields.io/pypi/v/dislocs?color=blue&logo=pypi&label=PyPI)](https://pypi.org/project/dislocs/)
[![Downloads](https://img.shields.io/pypi/dm/dislocs?color=green&logo=python&label=Downloads)](https://pypi.org/project/dislocs/)
[![License](https://img.shields.io/github/license/ZelongGuo/dislocs?logo=apache)](https://github.com/ZelongGuo/dislocs/blob/main/LICENSE)

# Disloc(ation)s

Calculating deformation, stress and strain from rectangular dislocation elements (RDE) and triangular
dislocation elements (TDE) sources in elastic half space.  
The implementations of this package is based on references (see References).


Mehdi's algorithm is recommended because of numerical instabilities and no artefact singularities, even though I also included Meade's implementation in the package.  
Mind the unit: recommend the International System of Units. Though strain is a dimensionless quantity, specifying dislocation unit same with fault parameters would help avoid further transformation.
## Requirments:  
- Numpy

## Features:
- Unified the strike-slip, dip-slip, and tensile coordinate systems for special faults (vertical faults, horizontal faults)
- Integration of TDE and RDE implementations, vertexes of TDE are order-free, more user-friendly

![An example of RDE and TDEs output](./tests/displacement_field.png)

## Installation
### By Pypi
The easiest way of building and installing the package is by PYPI (Installation under a python virtual env is recommended):  
```python
pip install dislocs
import dislocs as dl
```

### Build and install from the source codes
Or using CMake to make the installation from the source codes directly. In the projrct directory, and also strongly recommend executing under a python virtual env, the package (the shared library `dislocs.so` file ) is supposed to be installed in your current Python env automatically:  
```sh
mkdir -p build
cd build
cmake -DBUILD_LOCAL=ON ..
make install
```

## References:  
- Okada, Y., 1992, Internal deformation due to shear and tensile faults in a half-space, Bull. Seism. Soc. Am., 82, 1018-1040.

- Brendan J. Meade, 2007, Algorithms for the calculation of exact displacements, strains, and stresses for triangular dislocation elements in a uniform elastic half space, Computer & Geosciences, 33, 1064-1075.

- Nikkhoo, M., Walter, T. R., 2015, Triangular dislocation: an analytical, artefact-free solution.  Geophysical Journal International, 201, 2, 1119-1141.

> :copyright: Zelong Guo, Potsdam  
zelong.guo@outlook.com



