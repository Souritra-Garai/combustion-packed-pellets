#ifndef __KINETICS__
#define __KINETICS__

#include <math.h>
#include <vector>

#define R 8.314

#define A       1.76123E-2
#define E_a     26.0E3
#define Q_r     118.4E3
#define MW_Prod 85.674939E-3

// Calculates rate constant using Arrhenius rate law
// T - Temperature in K
float calc_k(float const &T);

// Calculates the derivative of rate constant with respect to T
// T - Temperature in K
float calc_k_prime(float const &T);

// Updates value eta according to the implicit scheme (Ref. Section 4.2)
// eta - holds older value of eta
// T - Temperature in K
// Delta_t - size of time step
// eta is updated in place
void update_eta(float &eta, float const &T, float Delta_t);

// vectorized version of update_eta
void update_eta(std::vector<float> &eta, std::vector<float> &T, float Delta_t);

#endif