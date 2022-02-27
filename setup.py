from distutils.command.build_ext import build_ext
from glob import glob
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension
import os


__version__ = "0.0.1"
extra_args = []

#Compiler args
if os.name != 'nt':
    extra_args.append('-std=c++11')                 #compile with c++11 standard
    extra_args.append('-static')                    #Statically link

SRCFOLDERS = ["src","PYsrc"]
SRCFILES=[]
for folder in SRCFOLDERS:
    SRCFILES=SRCFILES+glob(folder+"/*.cpp")
SRCFILES = sorted(SRCFILES)

ext_modules = [
    Pybind11Extension(
        "APC",
        SRCFILES,                  # Sort source files for reproducibility
        include_dirs=["include","extern/include"],
        extra_compile_args=extra_args,
        extra_objects=["extern/lib/cspice.a"]     #statically link cspice library https://stackoverflow.com/questions/4597228/how-to-statically-link-a-library-when-compiling-a-python-module-extension
    ),
]

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