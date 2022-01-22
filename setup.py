from distutils.command.build_ext import build_ext
from glob import glob
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension
__version__ = "0.0.1"
ext_modules = [
    Pybind11Extension(
        "APC",
        sorted(glob("src/*.cpp")),  # Sort source files for reproducibility
        include_dirs=["include"],
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