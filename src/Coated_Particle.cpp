#include "Coated_Particle.hpp"

Coated_Particle::Coated_Particle(Substance C, Substance S, long double d, long double D)
{
    Core    = &C;
    Shell   = &S;

    Core_Diameter = d;
    Overall_Diameter = D;
}

long double Coated_Particle::Get_Density()
{
    long double d3 = pow(Core_Diameter, 3);
    long double D3 = pow(Overall_Diameter, 3);

    return (d3 * Core->Get_Density() + (D3 - d3) * Shell->Get_Density()) / D3;
}

long double Coated_Particle::Get_Heat_Capacity()
{
    long double d3_rho_c = pow(Core_Diameter, 3) * Core->Get_Density();
    long double D3_d3_rho_s = (pow(Overall_Diameter, 3) - pow(Core_Diameter, 3)) * Shell->Get_Density();

    return (d3_rho_c * Core->Get_Heat_Capacity() + D3_d3_rho_s * Shell->Get_Heat_Capacity()) / (d3_rho_c + D3_d3_rho_s);
}

long double Coated_Particle::Get_Heat_Conductivity()
{
    long double d3_rho_c = pow(Core_Diameter, 3) * Core->Get_Density();
    long double D3_d3_rho_s = (pow(Overall_Diameter, 3) - pow(Core_Diameter, 3)) * Shell->Get_Density();

    return (d3_rho_c * Core->Get_Heat_Conductivity() + D3_d3_rho_s * Shell->Get_Heat_Conductivity()) / (d3_rho_c + D3_d3_rho_s);
}

