#include "Kinetics.hpp"

std::pair<float, float> calc_reaction_terms(float eta, float omega, float T, float Delta_t)
{
    float exp_term = exp(E_a / (R * T));
    
    float numerator     = - DH_r * A * (2 - 2 * eta - omega * Delta_t);
    float denominator   = Delta_t * A + 2 * exp_term;
    
    float first_order_coeff = (2 * E_a * exp_term / (R * T*T)) * numerator / (denominator*denominator);

    return std::pair<float, float> ((numerator/denominator - T*first_order_coeff), first_order_coeff);
}

float omega_update(float w_prev, float T, float n, float Delta_t)
{
    return A * (2 - 2*n - Delta_t * w_prev) / (Delta_t * A + 2 * exp(E_a / (R * T)));
}

float eta_update(float n_prev, float w, float Delta_t)
{
    float retval = n_prev + w * Delta_t;
    
    return retval > 1 ? 1 : retval; 
}

void reaction_update(
    std::vector<float>::iterator omega,
    std::vector<float>::iterator temperature,
    std::vector<float>::iterator eta,
    float Delta_t,
    std::vector<float>::iterator omega_arr_end
    )
{
    for(; omega < omega_arr_end; omega++, temperature++, eta++)
    {
        *omega = omega_update(*omega, *temperature, *eta, Delta_t);
        
        *eta = eta_update(*eta, *omega, Delta_t);
    }
}