#ifndef __KINETICS__
#define __KINETICS__

#include <math.h>
#include <vector>

#define R 8.314

#define A       465l
#define E_a     34.7E3l
#define DH_r    -118.4E3l

std::pair<long double, long double> Calc_Reaction_Terms(long double Eta, long double Temperature, long double Delta_t);

long double Omega_Update(long double Eta, long double Temperature, long double Delta_t);

long double Eta_Update(long double Eta, long double Omega, long double Delta_t);

long double Calc_Omega(long double conversion, long double temperature);

#endif