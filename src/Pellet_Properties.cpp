#include "Thermo_Physical_Properties/Pellet_Properties.hpp"

Pellet_Properties::Pellet_Properties(Coated_Particle P, Substance F, long double phi)
{
    Density             = Calc_Density(phi, P.Get_Density(), F.Get_Density());
    Heat_Capacity       = Calc_Heat_Capacity(phi, P.Get_Density(), P.Get_Heat_Capacity(), F.Get_Density(), F.Get_Heat_Capacity());
    Heat_Conductivity   = Calc_Heat_Conductivity_Bruggemann_Model(phi, P.Get_Heat_Conductivity(), F.Get_Heat_Conductivity());
}

long double Pellet_Properties::Get_Density()
{
    return Density;
}

long double Pellet_Properties::Get_Heat_Capacity()
{
    return Heat_Capacity;
}

long double Pellet_Properties::Get_Heat_Conductivity()
{
    return Heat_Conductivity;
}

Reaction_Zone_Pellet_Properties::Reaction_Zone_Pellet_Properties(Coated_Particle R, Coated_Particle P, Substance F, long double phi)
{
    long double rho_P       = 0.5 * (R.Get_Density() + P.Get_Density());
    long double C_P         = 0.5 * (R.Get_Heat_Capacity() + P.Get_Heat_Capacity());
    long double lambda_P    = 0.5 * (R.Get_Heat_Conductivity() + P.Get_Heat_Conductivity());

    Density             = Calc_Density(phi, rho_P, F.Get_Density());
    Heat_Capacity       = Calc_Heat_Capacity(phi, rho_P, C_P, F.Get_Density(), F.Get_Heat_Capacity());
    Heat_Conductivity   = Calc_Heat_Conductivity_Bruggemann_Model(phi, lambda_P, F.Get_Heat_Conductivity());
}

long double Reaction_Zone_Pellet_Properties::Get_Density()
{
    return Density;
}

long double Reaction_Zone_Pellet_Properties::Get_Heat_Capacity()
{
    return Heat_Capacity;
}

long double Reaction_Zone_Pellet_Properties::Get_Heat_Conductivity()
{
    return Heat_Conductivity;
}

long double Calc_Density(long double phi_p, long double rho_p, long double rho_f)
{
    return phi_p * rho_p + (1 - phi_p) * rho_f;
}

long double Calc_Heat_Capacity(long double phi_p, long double rho_p, long double c_p, long double rho_f, long double c_f)
{
    long double phi_p_rho_p = phi_p * rho_p;
    long double phi_f_rho_f = (1 - phi_p) * rho_f;

    return (phi_p_rho_p * c_p + phi_f_rho_f * c_f) / (phi_p_rho_p + phi_f_rho_f);
}

long double Calc_Heat_Conductivity_Bruggemann_Model(long double phi_p, long double lambda_p, long double lambda_f)
{
    long double D = pow((lambda_p / lambda_f) * (3*phi_p - 1), 2) + pow(2-3*phi_p, 2) + 2 * (lambda_p / lambda_f) * (2 + 9*phi_p - 9*phi_p*phi_p);

    return (lambda_f/4) * ((3*phi_p - 1)*(lambda_p / lambda_f) + 2 - 3*phi_p + sqrt(D));
}

