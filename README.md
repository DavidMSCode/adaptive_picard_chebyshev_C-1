# adaptive_picard_chebyshev pythonbind11
This repository contains the python wrapper for the adaptive Picard-Chebyshev C++ code for propagating orbits around the Earth, using cartesian coordinates.

## Disclaimer
This tool is experimental and should not be run on any production environment. There are currently few safeguards in the C++ code to catch errors. Bad input can and will cause a segmentation fault in the C++ binary and it will take the python kernel with it. To prevent loss of data it's recommended to run this program only within its own Python virtual environment.

## Build Requirements
```
All operating systems require Python 3.6+ and Pip 10+

[Windows]  
Visual Studio Code 2017+

[Linux-Debian]  
python3-dev  
build-essential

[Linux-RedHat]  
python3-devel  
build-essential

[Mac OS]  
Xcode Command Line Tools
```

## Compiling
### Mac OS

1. Make sure XCode Command Line Tools are installed. run "xcode-select --install" if not.
2. clone this repository
3. cd into the repository
4. run "git checkout pybind11wrapper" to move to this branch
5. run "pip install ."
6. in python "import APCexample"
7. use test.ipynb in bin/ to test that the module is working (The "matrices/" folder must be in the same working directory)

### Linux

### Windows
