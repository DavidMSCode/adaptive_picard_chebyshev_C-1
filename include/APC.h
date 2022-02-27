/*
*  AUTHORS:          David Stanley (davidms4@illinois.edu)
*  DATE WRITTEN:     Feb 2022
*  LAST MODIFIED:    Feb 2022
*  AFFILIATION:      Department of Aerospace Engineering, University of Illincois Champaign-Urbana
*  DESCRIPTION:      Methods that are acessible from python and binding code
*  REFERENCE:        Woollands, R., and Junkins, J., "Nonlinear Differential Equation Solvers
*                    via Adaptive Picard-Chebyshev Iteration: Applications in Astrodynamics", JGCD, 2016.
*/

#ifndef __APCMAIN__
#define __APCMAIN__

#include <vector>
#include <Orbit.h>

std::vector<std::vector<double> > PropagateICs(std::vector<double> r, std::vector<double> v, double t0, double tf);

class Orbit PropagateOrbit(std::vector<double> r, std::vector<double> v, double t0, double tf);

#endif