#include "Kinetics.hpp"

std::pair<long double, long double> Calc_Reaction_Terms(long double eta, long double T, long double Delta_t)
{
    // long double exp_term = exp(E_a / (R * T));
    
    // long double numerator     = - DH_r * A * (1.0l - eta);
    // long double denominator   = A * Delta_t + exp_term;
    
    // long double first_order_coeff = (E_a * exp_term / (R * T*T)) * numerator / (denominator*denominator);

    // return std::pair<long double, long double> ((numerator/denominator - T*first_order_coeff), first_order_coeff);

    long double omega = - DH_r * A*exp(-E_a/(R*T));

    long double first_order_coeff = omega * E_a / (R * T * T);

    return std::pair<long double, long double> (omega - T * first_order_coeff, first_order_coeff);
}

long double Omega_Update(long double eta, long double T, long double Delta_t)
{
    return A * (1 - eta) / (A * Delta_t + exp(E_a / (R * T)));
}

long double Eta_Update(long double eta, long double omega, long double Delta_t)
{
    long double retval = eta + omega * Delta_t;
    
    return std::max(0.0l, std::min(1.0l, retval)); 
}

long double Calc_Omega(long double eta, long double T)
{
    return A * exp(- E_a / (R*T)) * (1.0l - eta);
}