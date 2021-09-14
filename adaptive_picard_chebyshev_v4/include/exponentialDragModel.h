
#ifndef DRAG_MODEL_H
#define DRAG_MODEL_H
#include "satellite_properties.h"

double atmospheric_density(double alt);
int indexSearch(double *p, int length_t, double key);
void drag_acceleration(double t, double* X, double* V, double* aECEF,struct satellite_properties sat);

#endif