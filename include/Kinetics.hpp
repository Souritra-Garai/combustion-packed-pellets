#ifndef __KINETICS__
#define __KINETICS__

#include <math.h>

#define c 6
#define r 32.5E-6

#define D_0 9.54E-8
#define E_A 26E3

#define R 8.314

float get_D(float T);

float get_t_b(float T);

float get_omega(float T, float eta);

#endif