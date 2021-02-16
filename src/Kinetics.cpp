#include "Kinetics.hpp"

float calc_k(float const &T)
{
    return A * exp(- E_a / (R*T));
}

void update_eta(float &eta, float const &T, float Delta_t)
{
    float k = calc_k(T);
    eta = (eta + k * Delta_t) / (1 + k * Delta_t);

    if (eta > 1) eta = 1;
    if (eta < 0) eta = 0;
}

void update_eta(std::vector<float> &eta, std::vector<float> &T, float Delta_t)
{
    if (eta.size() != T.size()) throw "Input array sizes eta and T are dissimilar!!";

    for (std::vector<float>::iterator eta_i = eta.begin(), T_i = T.begin(); eta_i < eta.end(); eta_i++, T_i++)

        update_eta(*eta_i, *T_i, Delta_t);
}