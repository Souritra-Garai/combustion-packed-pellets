#include "Kinetics.hpp"

float get_D(float T)
{
    return D_0 * exp(- E_A / (R * T));
}

float get_t_b(float T)
{
    return pow(r, 2) / (c * get_D(T));
}

float get_omega(float T, float eta)
{
    if (eta > 1) throw "eta is greater than 1";

    if (eta <= 0.0) return 0.0;

    return (1 - eta) / get_t_b(T);    
}