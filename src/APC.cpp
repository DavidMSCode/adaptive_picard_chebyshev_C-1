/*
*  AUTHORS:          David Stanley (davidms4@illinois.edu)
*  DATE WRITTEN:     Feb 2022
*  LAST MODIFIED:    Feb 2022
*  AFFILIATION:      Department of Aerospace Engineering, University of Illincois Champaign-Urbana
*  DESCRIPTION:      Methods that are acessible from python and binding code
*  REFERENCE:        Woollands, R., and Junkins, J., "Nonlinear Differential Equation Solvers
*                    via Adaptive Picard-Chebyshev Iteration: Applications in Astrodynamics", JGCD, 2016.
*/

#include <adaptive_picard_chebyshev.h>
#include <c_functions.h>
#include <EGM2008.h>
#include <time.h> 
#include <errno.h>
#include <vector>
#include <Orbit.h>
#include <APC.h>

std::vector<std::vector<double> > PropagateICs(std::vector<double> r, std::vector<double> v, double t0, double tf){
  printf("%s",typeid(r).name());
  //Convert vectors to array since pybind wants vectors but the functions are coded for arrays
  double* r0 = &r[0];
  double* v0 = &v[0];
  double dt    = 30.0;                             // Soution Output Time Interval (s)
  double deg   = 70.0;                             // Gravity Degree (max 100)
  double tol   = 1.0e-15;                          // Tolerance

  // Initialize Output Variables
  int soln_size = int(1.1*(tf/dt));
  if (soln_size == 1){
    soln_size = 2;
  }
  std::vector<double> Soln(soln_size*6,0.0);
  //Soln = static_cast<double*>(calloc(soln_size*6,sizeof(double)));       // Position (km) & Velocity (km/s)

  double Feval[2] = {0.0};

  // Call Adaptive Picard Chebyshev Integrator
  clock_t startTime = clock();
  std::vector<std::vector<double> > states;
  // for (int tt=0; tt<=1; tt++){
  states = adaptive_picard_chebyshev(r0,v0,t0,tf,dt,deg,tol,soln_size,Feval,Soln);
  // }
  clock_t endTime = clock();
  float elapsedTime = ((float) (endTime - startTime))/CLOCKS_PER_SEC/(1.0);
  printf("Elapsed time: %f s\t",elapsedTime);

  // Number of function evaluations
  int total;
  total = int(ceil(Feval[0] + Feval[1]*pow(6.0,2)/pow(deg,2)));
  printf("Func Evals: %i\t",total);

  // Save as text file [time r v H]
  std::vector<std::vector<double> > States(6);
  double state[6] = {0.0};
  std::vector<double> Hs;
  std::vector<double> Ts;
  double H    = 0.0;
  double H0   = 0.0;
  double Hmax = 0.0;

  double t_curr = t0;
  for (int i=1; i<=soln_size; i++){
    Ts.push_back(t_curr);
    for (int j=1; j<=6; j++){
      States[j-1].push_back(Soln[ID2(i,j,soln_size)]);
      state[j-1] = Soln[ID2(i,j,soln_size)];
    }
    jacobiIntegral(t_curr,state,&H,deg);
    if (i == 1){
      H0 = H;
    }
    if (fabs((H-H0)/H0) > Hmax){
      Hmax = fabs((H-H0)/H0);
    }
    Hs.push_back(fabs((H-H0)/H0));
    t_curr = t_curr + dt;
    if (t_curr > tf){
      break;
    }
  }
  printf("Hmax %1.16E\n",Hmax);
  //Assemble solution vector
  std::vector<std::vector<double> > Solution;
  Solution.push_back(Ts);
  for(int i=0; i<=5; i++){
    Solution.push_back(States[i]);
  }
  Solution.push_back(Hs);
  //free(Soln);
  return Solution;
}

class Orbit PropagateOrbit(std::vector<double> r, std::vector<double> v, double t0, double tf){
  std::vector<std::vector<double> > solution;
  solution = PropagateICs(r,v,t0,tf);
  Orbit orbit(solution);
  return orbit;
}

