#ifndef __KINETICS__
#define __KINETICS__

#include <math.h>
#include <vector>

#define R 8.314

#define A       0 // 2.215E25
#define E_a     448.4E3
#define DH_r    -50E3

std::pair<long double, long double> calc_reaction_terms(long double eta, long double omega, long double T, long double Delta_t);

long double omega_update(long double omega_prev, long double temperature, long double eta, long double Delta_t);

long double eta_update(long double omega, long double eta_prev, long double Delta_t);

void reaction_update(
    std::vector<long double>::iterator omega_prev,
    std::vector<long double>::iterator temperature,
    std::vector<long double>::iterator eta,
    long double Delta_t,
    std::vector<long double>::iterator omega_arr_end
    );

long double calc_conversion_rate(long double conversion, long double temperature);

#endif