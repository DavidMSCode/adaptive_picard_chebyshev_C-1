from distutils.command.build_ext import build_ext
from glob import glob
import platform
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension
import os


__version__ = "0.0.1"
cpp_extra_args = []

"""Compiler args"""
if os.name != 'nt':
    cpp_extra_args.append('-std=c++11')                             #compile with c++11 standard if not Windows
# if platform.system() == 'Darwin':
#    cpp_extra_args.append('-mmacosx-version-min=10.9')

"""Get all necessary C++ source files"""
SRCFOLDERS = ["src","PYsrc"]
SRCFILES=[]
for folder in SRCFOLDERS:
    SRCFILES=SRCFILES+glob(folder+"/*.cpp")
SRCFILES = sorted(SRCFILES)                                         # Sort source files for consistency

"""Define external module"""
ext_modules = [
    Pybind11Extension(
        "APC",
        SRCFILES,                                                   #.cpp files to be compiled
        include_dirs=["include","extern/cspice/include"],           #header file locations for src files and external libraries
        extra_compile_args=cpp_extra_args,                          #extra cpp compile args
        extra_objects=["extern/cspice/lib/cspice.a"]                #statically link cspice library https://stackoverflow.com/questions/4597228/how-to-statically-link-a-library-when-compiling-a-python-module-extension
    ),
]

"""Setup module"""
setup(
    name = "APC",
    version = __version__,
    cmdclass={"build_ext": build_ext},
    author="David Stanley",
    author_email = "davidms4@illinois.edu",
    description = "Test compilation",
    ext_modules = ext_modules,
    zip_safe = False,
    python_requires = ">=3.6",
)