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

    Core_Diameter = d;
    Overall_Diameter = D;

    Moles_of_Core_Material = (M_PIl / 6) * d3_rho_c / C.Get_Molecular_Weight();
    Moles_of_Shell_Material = (M_PIl / 6) * D3_d3_rho_s / S.Get_Molecular_Weight();
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

void Coated_Particle::Write_to_File(std::ofstream &file, const char *name)
{
    file << name << " Properties" << std::endl;

    file << std::endl;

    file << "Core Diameter:\t" << Core_Diameter << "\t m" << std::endl;
    file << "Overall Diameter:\t" << Overall_Diameter << "\t m" << std::endl;

    file << std::endl;

    file << "Moles of Core Material:\t" << Moles_of_Core_Material << "\t mol." << std::endl;
    file << "Moles of Shell Material:\t" << Moles_of_Shell_Material << "\t mol." << std::endl;

    file << std::endl;

    file << "Density :\t" << Density << "\t kg / m3" << std::endl;
    file << "Heat Capacity :\t" << Heat_Capacity << "\t J / kg - K" << std::endl;
    file << "Heat Conductivity :\t" << Heat_Conductivity << "\t W / m - K" << std::endl;

    file << std::endl;
}

long double Coated_Particle::Get_Moles_of_Core_Material()
{
    return Moles_of_Core_Material;
}

long double Coated_Particle::Get_Moles_of_Shell_Material()
{
    return Moles_of_Shell_Material;
}

long double Coated_Particle::Get_Core_Material_Concentration()
{
    return Moles_of_Core_Material / (M_PIl * pow(Overall_Diameter, 3) / 6);
}

long double Coated_Particle::Get_Shell_Material_Concentration()
{
    return Moles_of_Shell_Material / (M_PIl * pow(Overall_Diameter, 3) / 6);
}

long double Calc_Product_Particle_Core_Diameter(long double n, Substance s)
{
    return pow( n * s.Get_Molecular_Weight() / (s.Get_Density() * (M_PIl / 6)), 1.0/3.0 );
}

long double Calc_Product_Particle_Overall_Diameter(long double n, Substance s, long double d)
{
    return pow( pow(d, 3) + n * s.Get_Molecular_Weight() / (s.Get_Density() * (M_PIl / 6)), 1.0/3.0 );
}
