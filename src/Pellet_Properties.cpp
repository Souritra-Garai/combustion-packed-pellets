#include "Pellet_Properties.hpp"

Pellet_Properties::Pellet_Properties(Coated_Particle P, Substance F, long double phi)
{
    Particle_A = &P;
    Fluid = &F;

    Packing_Volume_Fraction = phi;

    Density             = Calc_Density(Packing_Volume_Fraction, Particle_A->Get_Density(), Fluid->Get_Density());
    Heat_Capacity       = Calc_Heat_Capacity(Packing_Volume_Fraction, Particle_A->Get_Density(), Particle_A->Get_Heat_Capacity(), Fluid->Get_Density(), Fluid->Get_Heat_Capacity());
    Heat_Conductivity   = Calc_Heat_Conductivity_Bruggemann_Model(Packing_Volume_Fraction, Particle_A->Get_Heat_Conductivity(), Fluid->Get_Heat_Conductivity());
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

Reaction_Zone_Pellet_Properties::Reaction_Zone_Pellet_Properties(Coated_Particle R, Coated_Particle P, Substance F, long double phi) : Pellet_Properties(R, F, phi)
{
    Particle_B = &P;

    Density             = Calc_Density(Packing_Volume_Fraction, 0.5 * (Particle_A->Get_Density() + Particle_B->Get_Density()), Fluid->Get_Density());
    Heat_Capacity       = Calc_Heat_Capacity(Packing_Volume_Fraction, 0.5 * (Particle_A->Get_Density() + Particle_B->Get_Density()), 0.5* (Particle_A->Get_Heat_Capacity() + Particle_B->Get_Heat_Capacity()), Fluid->Get_Density(), Fluid->Get_Heat_Capacity());
    Heat_Conductivity   = Calc_Heat_Conductivity_Bruggemann_Model(Packing_Volume_Fraction, 0.5*(Particle_A->Get_Heat_Conductivity()+Particle_B->Get_Heat_Conductivity()), Fluid->Get_Heat_Conductivity());
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

