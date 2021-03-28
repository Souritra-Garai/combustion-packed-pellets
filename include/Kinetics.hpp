#ifndef __KINETICS__
#define __KINETICS__

#include <math.h>
#include <vector>

#define R 8.314

#define A       2.215E25
#define E_a     448.4E3
#define DH_r    -50E3

std::pair<float, float> reaction_term(float eta, float omega, float T, float Delta_t);

float omega_update(float omega_prev, float temperature, float eta, float Delta_t);

float eta_update(float omega, float eta_prev, float Delta_t);

void reaction_update(
    std::vector<float>::iterator omega_prev,
    std::vector<float>::iterator temperature,
    std::vector<float>::iterator eta,
    float Delta_t,
    std::vector<float>::iterator omega_arr_end
    );

#endif