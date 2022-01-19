from glob import glob
import sys
from pybind11 import get_cmake_dir
# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

__version__ = "0.0.1"
ext_modules = [
    Pybind11Extension(
        "APCexample",
        sorted(glob("./src/*.cpp")),  # Sort source files for reproducibility
        include_dirs=["./include"],
        # Example: passing in the version to the compiled code
        #efine_macros = [('VERSION_INFO', __version__)],
    ),
]

setup(
    name = "APCexample",
    version = __version__,
    author="David Stanley",
    author_email = "davidms4@illinois.edu",
    description = "Test compilation",
    long_description="",
    ext_modules=ext_modules,
    # Currently, build_ext only provides an optional "highest supported C++
    # level" feature, but in the future it may provide more features.
    zip_safe=False,
    python_requires=">=3.6",
)