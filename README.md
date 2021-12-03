# adaptive_picard_chebyshev
This repository contains the python wrapper for the adaptive Picard-Chebyshev C++ code for propagating orbits around the Earth, using cartesian coordinates.

To install this module on Mac OS

1. Make sure XCode Command Line Tools are installed. run "xcode-select --install" if not.
2. clone this repository
3. cd into the repository
4. run "git checkout pybind11wrapper" to move to this branch
5. run "pip install ."
6. in python "import APCexample"
7. use test.ipynb in bin/ to test that the module is working (The "matrices/" folder must be in the same working directory)
