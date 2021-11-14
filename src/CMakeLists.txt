cmake_minimum_required(VERSION 3.0.0)
project(APyC VERSION 0.1.0)

#Find the most up to date python package
find_package(Python COMPONENTS Interpreter Development)

#Header files
include_directories(../include)
#Includes for Pybind11
include_directories(extern/pybind11/include/pybind11)
add_subdirectory(extern/pybind11)

#Compile module using pybind11 function
pybind11_add_module(APyC APyC.cpp adaptive_picard_chebyshev.cpp polydegree_segments.cpp picard_chebyshev_propagator.cpp prepare_propagator.cpp
picard_iteration.cpp clenshaw_curtis_ivpII.cpp rv2elm.cpp c_functions.cpp FandG.cpp perigee_approx.cpp EGM2008.cpp
lsq_chebyshev_fit.cpp chebyshev.cpp radial_gravity.cpp matrix_loader.cpp eci2ecef.cpp ecef2eci.cpp perturbed_gravity.cpp
picard_error_feedback.cpp interpolate.cpp reosc_perigee.cpp)