/*
*  AUTHORS:          David Stanley (davidms4@illinois.edu)
*  DATE WRITTEN:     Feb 2022
*  LAST MODIFIED:    Feb 2022
*  AFFILIATION:      Department of Aerospace Engineering, University of Illincois Champaign-Urbana
*  DESCRIPTION:      Class that stores orbit solution and properties
*  REFERENCE:        Woollands, R., and Junkins, J., "Nonlinear Differential Equation Solvers
*                    via Adaptive Picard-Chebyshev Iteration: Applications in Astrodynamics", JGCD, 2016.
*
*/


#include "Orbit.h"
#include <vector>


//Orbit Constructors
Orbit::Orbit(std::vector<std::vector<double > > Solution){
    SetSolution(Solution);
}
Orbit::Orbit(double area, double reflectivity, double mass){
    SetProperties(area, reflectivity, mass);
}

//Setter Functions
void Orbit::SetSolution(std::vector<std::vector<double > > Solution){
    Soln = Solution;
}

void Orbit::SetProperties(double area, double reflectivity, double mass){
    properties.Area = area;
    properties.Mass = mass;
    properties.Reflectivity = reflectivity;
}
