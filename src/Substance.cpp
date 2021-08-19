#include "Thermo_Physical_Properties/Substance.hpp"

Substance::Substance(long double rho, long double C, long double lambda, long double MW)
{
    Density = rho;
    Heat_Capacity = C;
    Heat_Conductivity = lambda;

    Molecular_Weight = MW;
}

long double Substance::Get_Density()            {return Density;}
long double Substance::Get_Heat_Capacity()      {return Heat_Capacity;}
long double Substance::Get_Heat_Conductivity()  {return Heat_Conductivity;}

long double Substance::Get_Molecular_Weight()  {return Molecular_Weight;}

void Substance::Write_to_File(std::ofstream &file, const char *name)
{
    file << name << " Properties" << std::endl;

    file << "Density :\t" << Get_Density() << "\t kg / m3" << std::endl;
    file << "Heat Capacity :\t" << Get_Heat_Capacity() << "\t J / kg - K" << std::endl;
    file << "Heat Conductivity :\t" << Get_Heat_Conductivity() << "\t W / m - K" << std::endl;

    file << std::endl;
}