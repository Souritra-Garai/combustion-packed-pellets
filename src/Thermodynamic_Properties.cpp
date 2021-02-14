#include "Thermodynamic_Properties.hpp"

float get_lambda_m(float lambda_p, float lambda_f, float phi_p, float alpha_p_1)
{
    if (phi_p > 1 || phi_p < 0) throw "phi_p is not within [0,1]!!";
    
    if (alpha_p_1 > 1 || alpha_p_1 < 0) throw "alpha_p_1 is not within [0,1]!!";

    float phi_f = 1 - phi_p;

    float D = (2*lambda_p - lambda_f) * phi_p * (1 - alpha_p_1) + (2*lambda_f - lambda_p) * ( (phi_f + phi_p*alpha_p_1 - 0.5) / phi_f );

    return 0.5 * ( D + sqrt(pow(D, 2) + 2*lambda_p*lambda_f) );
}

float get_rho_m(float phi_p, float rho_f)
{
    if (phi_p > 1 || phi_p < 0) throw "phi_p is not within [0,1]!!";

    return phi_p * rho_p + (1-phi_p) * rho_f;
}

float get_Cp_m(float phi_p, float Cp_p, float Cp_f, float rho_f)
{
    if (phi_p > 1 || phi_p < 0) throw "phi_p is not within [0,1]!!";
    
    float X_p = phi_p * rho_p / ( phi_p * rho_p + (1 - phi_p) * rho_f );

    return X_p * Cp_p + (1 - X_p) * Cp_f;
}