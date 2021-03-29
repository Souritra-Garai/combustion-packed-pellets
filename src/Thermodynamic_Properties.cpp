#include "Thermodynamic_Properties.hpp"

long double calc_lambda_m(long double phi_p, long double alpha_p_ME, long double lambda_p, long double lambda_f)
{
    if (phi_p > 1 || phi_p < 0) throw "phi_p is not within [0,1]!!";
    
    if (alpha_p_ME > 1 || alpha_p_ME < 0) throw "alpha_p_1 is not within [0,1]!!";

    long double phi_f = 1 - phi_p;

    long double D = (2*lambda_p - lambda_f) * phi_p * (1 - alpha_p_ME) + (2*lambda_f - lambda_p) * ( (phi_f + phi_p*alpha_p_ME - 0.5) / phi_f );

    return 0.5 * ( D + sqrt(pow(D, 2) + 2*lambda_p*lambda_f) );
}

long double calc_rho_m(long double phi_p, long double rho_p, long double rho_f)
{
    if (phi_p > 1 || phi_p < 0) throw "phi_p is not within [0,1]!!";

    return phi_p * rho_p + (1 - phi_p) * rho_f;
}

long double calc_Cp_m(long double phi_p, long double rho_p, long double rho_f, long double Cp_p, long double Cp_f)
{
    if (phi_p > 1 || phi_p < 0) throw "phi_p is not within [0,1]!!";
    
    return ( phi_p * rho_p * Cp_p + (1 - phi_p) * rho_f * Cp_f ) / ( phi_p * rho_p + (1 - phi_p) * rho_f );
}