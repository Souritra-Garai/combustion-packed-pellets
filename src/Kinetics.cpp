#include "Kinetics.hpp"

std::pair<long double, long double> calc_reaction_terms(long double eta, long double omega, long double T, long double Delta_t)
{
    long double exp_term = exp(E_a / (R * T));
    
    long double numerator     = - DH_r * A * (2 - 2 * eta - omega * Delta_t);
    long double denominator   = Delta_t * A + 2 * exp_term;
    
    long double first_order_coeff = (2 * E_a * exp_term / (R * T*T)) * numerator / (denominator*denominator);

    return std::pair<long double, long double> ((numerator/denominator - T*first_order_coeff), first_order_coeff);
}

long double omega_update(long double w_prev, long double T, long double n, long double Delta_t)
{
    return A * (2 - 2*n - Delta_t * w_prev) / (Delta_t * A + 2 * exp(E_a / (R * T)));
}

long double eta_update(long double n_prev, long double w, long double Delta_t)
{
    long double retval = n_prev + w * Delta_t;
    
    return std::max(0.0l, std::min(1.0l, retval)); 
}

void reaction_update(
    std::vector<long double>::iterator omega,
    std::vector<long double>::iterator temperature,
    std::vector<long double>::iterator eta,
    long double Delta_t,
    std::vector<long double>::iterator omega_arr_end
    )
{
    for(; omega < omega_arr_end; omega++, temperature++, eta++)
    {
        *omega = omega_update(*omega, *temperature, *eta, Delta_t);
        
        *eta = eta_update(*eta, *omega, Delta_t);
    }
}

long double calc_conversion_rate(long double eta, long double T)
{
    return A * exp(- E_a / (R*T)) * (1 - std::max(0.0l, std::min(1.0l, eta)));
}