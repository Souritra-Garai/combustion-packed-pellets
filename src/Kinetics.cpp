#include "Kinetics.hpp"

std::pair<float, float> reaction_term(float eta, float omega, float T, float Delta_t)
{
    float exp_term = exp(E_a / (R * T));
    
    float numerator     = A * (2 - 2 * eta - omega * Delta_t);
    float denominator   = Delta_t * A + 2 * exp_term;
    
    float first_order_coeff = (2 * E_a * exp_term / (R * T*T)) * numerator / (denominator*denominator);

    return std::pair<float, float> ((numerator/denominator - T*first_order_coeff), first_order_coeff);
}