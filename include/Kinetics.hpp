#ifndef __KINETICS__
#define __KINETICS__

#include <math.h>
#include <vector>

#define R 8.314

#define A       2.215E25
#define E_a     448.4E3
#define Q_r     118.4E3
#define MW_Prod 85.674939E-3

std::pair<float, float> reaction_term(float eta, float omega, float T, float Delta_t);

#endif