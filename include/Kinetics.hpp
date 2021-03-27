#ifndef __KINETICS__
#define __KINETICS__

#include <math.h>
#include <vector>

#define R 8.314

#define A       2.215E25
#define E_a     448.4E3
#define DH_r    -50E3

std::pair<float, float> reaction_term(float eta, float omega, float T, float Delta_t);

#endif