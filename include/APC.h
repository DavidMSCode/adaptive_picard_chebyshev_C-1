/*
*  AUTHORS:          David Stanley (davidms4@illinois.edu)
*  DATE WRITTEN:     Feb 2022
*  LAST MODIFIED:    Feb 2022
*  AFFILIATION:      Department of Aerospace Engineering, University of Illinois Urbana-Champaign
*  DESCRIPTION:      Methods that are acessible from python and binding code
*  REFERENCE:        Woollands, R., and Junkins, J., "Nonlinear Differential Equation Solvers
*                    via Adaptive Picard-Chebyshev Iteration: Applications in Astrodynamics", JGCD, 2016.
*/

#ifndef __APCMAIN__
#define __APCMAIN__

#include <vector>
#include <Orbit.h>

std::vector<std::vector<double> > PropagateICs(std::vector<double> r, std::vector<double> v, double t0, double tf, double area, double reflectance, double mass, double drag_C, bool compute_drag, bool compute_srp, bool compute_third_body);

class Orbit PropagateOrbit(std::vector<double> r, std::vector<double> v, double t0, double tf, double area, double reflectance, double mass, double drag_C, bool compute_drag, bool compute_srp, bool compute_third_body);

#endif