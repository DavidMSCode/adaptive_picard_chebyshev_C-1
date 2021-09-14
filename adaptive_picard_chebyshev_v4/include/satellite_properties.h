//This header defines the satellite_properties struct which carries properties for calculating drag

#ifndef SATELLITE_PROPERTIES_H
#define SATELLITE_PROPERTIES_H

struct satellite_properties
{
    double A;       //Area of satellite (m^2)
    double cd;      //drag coefficient of satellite
    double m;       //mass of satellite (kg)
};

#endif