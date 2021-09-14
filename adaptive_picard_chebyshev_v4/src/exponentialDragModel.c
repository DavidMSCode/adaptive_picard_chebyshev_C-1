/*
*  AUTHORS:          David Stanley (davidms4@illinois.edu)
*  DATE WRITTEN:     Aug 2021
*  LAST MODIFIED:    Aug 2021
*  AFFILIATION:      Department of Aerospace Engineering, University of Illinois Urbana-Champaign, Champaign, IL
*  DESCRIPTION:      Atmospheric drag functions utilizing Vallado's exponential atmospheric model
*
* 
* REFERENCES:
* 1. Vallado, David A, and Wayne D McClain. Fundamentals of Astrodynamics and Applications. 3rd ed., Springer, 2007.
*
* COMMENTS:
*
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../include/exponentialDragModel.h"
#include "const.h"


double density_vector[28] = {
    1.225,
    3.899e-2,
    1.774e-2,
    3.972e-3,
    1.057e-3,
    3.206e-4,
    8.770e-5,
    1.905e-5,
    3.396e-6,
    5.297e-7,
    9.661e-8,
    2.438e-8,
    8.484e-9,
    3.845e-9,
    2.070e-9,
    5.464e-10,
    2.789e-10,
    7.248e-11,
    2.418e-11,
    9.518e-12,
    3.725e-12,
    1.585e-12,
    6.967e-13,
    1.454e-13,
    3.614e-14,
    1.170e-14,
    5.245e-15,
    3.019e-15};

    double base_alt_vector[28] = {
    0,
    25,
    30,
    40,
    50,
    60,
    70,
    80,
    90,
    100,
    110,
    120,
    130,
    140,
    150,
    180,
    200,
    250,
    300,
    350,
    400,
    450,
    500,
    600,
    700,
    800,
    900,
    1000};
  
    double scale_height_vector[28] = {
    7.249,
    6.349,
    6.682,
    7.554,
    8.382,
    7.714,
    6.549,
    5.799,
    5.382,
    5.877,
    7.263,
    9.473,
    12.636,
    16.149,
    22.523,
    29.740,
    37.105,
    45.546,
    53.628,
    53.298,
    58.515,
    60.828,
    63.822,
    71.835,
    88.667,
    124.64,
    181.05,
    268.00};

int indexSearch(double *p, int length_t, double key){
    /*This function attempts to find the index for the value given by the key by searching from largest
     to smallest. If the key is not present the function returns the index of the first item smaller 
     than the key value assuming a sorted list. A binary search would probably be more efficient, 
     but the arrays are very short.
    */
    //Set flags and index for iterating through vector
    bool indexFound=false;  
    int i = length_t-1;
    int index_sol=-1;

    while(!indexFound && i>=0){
        if (*p <= key){
            //printf("%f<=%f\n",*p,key);
            //Smaller value found. End loop.
            index_sol=i;
            indexFound=true;
        }
        //increment search index and array pointer
        i--;
        p--;
    }
    //Return found index. If value is -1 thenn an index was not found.
    return index_sol;
}

double atmospheric_density(double alt){
    //Find the index for the base altitude and get the corresponding base altitude, scaling height and atmospheric density
    int length_t = sizeof(base_alt_vector)/sizeof(base_alt_vector[0]);                   //Get length of array
    double *p = (double *)(&base_alt_vector + 1) - 1;                                    //Get pointer for last element in array
    int index = indexSearch(p, length_t, alt);                                           //Find base altitude index
    //printf("index: %d\n",index);

    double base_density = density_vector[index];
    double base_alt = base_alt_vector[index];
    double scale_height = scale_height_vector[index];
    //printf("rho_0:%E\nalt_0:%f\nH:%f\n",base_density,base_alt,scale_height);
    //calculate atmospheric density at given altitude
    double density = base_density*exp(-(alt-base_alt)/scale_height);
    return density;
}

void drag_acceleration(double t, double* X, double* V, double* dragECEF, struct satellite_properties sat){
    
    //Takes the satellite mass(kg), drag coefficient, area(m^2), position(km) and velocity(km/s) and calculates the acceleration (km/s^2) due to atmospheric drag.
    double r = sqrt(pow(X[0],2) + pow(X[1],2) + pow(X[2],2));    //distance from origin in ECEF km
    double v = sqrt(pow(V[0],2) + pow(V[1],2) + pow(V[2],2));    //speed in ECEF km/s
    double A = sat.A/pow(1000,2);   //Sat area km^2
    double cd = sat.cd;             //Sat drag coefficient
    double mass = sat.m;            //Sat mass kg
    double r_eq = C_Req;           //radius of the Earth (km)
    double alt = r - r_eq;          //altitude of satellite above the Earth (km)

    //Get the atmospheric density kg/km^3
    double rho = atmospheric_density(alt)*pow(1000,3);
    
    //Calculate drag acceleration km/s^2
    double drag_a = 0.5*rho*v*A*cd/mass;
    
    //Get acceleration vector in the direction opposite of the velocity vector
    for (int i=0;i<3;i++){
        dragECEF[i]= -drag_a*V[i];
    }
}