#include "Thermo_Physical_Properties/Coated_Particle.hpp"

Coated_Particle::Coated_Particle(Substance C, Substance S, long double d, long double D)
{
    long double d3 = pow(d, 3);
    long double D3 = pow(D, 3);

    Density = (d3 * C.Get_Density() + (D3 - d3) * S.Get_Density()) / D3;

    long double d3_rho_c = d3 * C.Get_Density();
    long double D3_d3_rho_s = (D3 - d3) * S.Get_Density();

    Heat_Capacity = (d3_rho_c * C.Get_Heat_Capacity() + D3_d3_rho_s * S.Get_Heat_Capacity()) / (d3_rho_c + D3_d3_rho_s);

    Heat_Conductivity = (d3_rho_c * C.Get_Heat_Conductivity() + D3_d3_rho_s * S.Get_Heat_Conductivity()) / (d3_rho_c + D3_d3_rho_s);
}

long double Coated_Particle::Get_Density()
{
    return Density;
}

long double Coated_Particle::Get_Heat_Capacity()
{
    return Heat_Capacity;
}

long double Coated_Particle::Get_Heat_Conductivity()
{
    return Heat_Conductivity;
}

