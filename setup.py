from glob import glob
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension
__version__ = "0.0.1"
cpp_args = ['-std=c++11']
ext_modules = [
    Pybind11Extension(
        "APCexample",
        sorted(glob("src/*.cpp")),  # Sort source files for reproducibility
        include_dirs=["include"],
        extra_compile_args = cpp_args,
        cxx_std = 11
    ),
]

setup(
    name = "APCexample",
    version = __version__,
    author="David Stanley",
    author_email = "davidms4@illinois.edu",
    description = "Test compilation",
    ext_modules = ext_modules,
    zip_safe = False,
    python_requires = ">=3.6",
)