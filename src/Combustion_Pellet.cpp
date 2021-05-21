#include "Thermo_Physical_Properties/Combustion_Pellet.hpp"

// #include <iostream>

Combustion_Pellet::Combustion_Pellet(
    long double D,
    long double L,
    Pellet_Properties PH,
    Pellet_Properties PC,
    Reaction_Zone_Pellet_Properties R
) : Pre_Heat_Zone(PH), Post_Combustion_Zone(PC), Reaction_Zone(R)
{
    Diameter = D;
    Length = L;

    Convective_Heat_Transfer_Coefficient = 0;
    Radiative_Emissivity = 0;
    Ambient_Temperature = 297;

    Ignition_Temperature = 300;

    // std::cout << "Density\t" << Pre_Heat_Zone.Get_Density() << '\t' << Reaction_Zone.Get_Density() << '\t' << Post_Combustion_Zone.Get_Density() << std::endl;
    // std::cout << "Heat Capacity\t" << Pre_Heat_Zone.Get_Heat_Capacity() << '\t' << Reaction_Zone.Get_Heat_Capacity() << '\t' << Post_Combustion_Zone.Get_Heat_Capacity() << std::endl;
    // std::cout << "Heat Conductivity\t" << Pre_Heat_Zone.Get_Heat_Conductivity() << '\t' << Reaction_Zone.Get_Heat_Conductivity() << '\t' << Post_Combustion_Zone.Get_Heat_Conductivity() << std::endl;
}

void Combustion_Pellet::Set_Convective_Heat_Transfer_Coefficient(long double H)
{
    Convective_Heat_Transfer_Coefficient = H;
}

void Combustion_Pellet::Set_Radiative_Emissivity(long double e)
{
    Radiative_Emissivity = e;
}

void Combustion_Pellet::Set_Ambient_Temperature(long double T)
{
    Ambient_Temperature = T;
}

void Combustion_Pellet::Set_Ignition_Temperature(long double T)
{
    Ignition_Temperature = T;
}

long double Combustion_Pellet::Get_Pellet_Length()
{
    return Length;
}

void Combustion_Pellet::Write_to_File(std::ofstream &file, const char *name)
{
    file << "Properties of " << name << std::endl;

    file << "Diameter :\t" << Diameter << "\t m" << std::endl;
    file << "Length :\t" << Length << "\t m" << std::endl;

    file << std::endl;

    Pre_Heat_Zone.Write_to_File(file, "Pre-Heat Zone Pellet");
    Reaction_Zone.Write_to_File(file, "Reaction Zone Pellet");
    Post_Combustion_Zone.Write_to_File(file, "Post-Combustion Zone Pellet");

    file << "Pellet Environment" << std::endl;
    file << "Ambient Temperature :\t" << Ambient_Temperature << "\t K" << std::endl;
    file << "Radiative Emissivity :\t" << Radiative_Emissivity << std::endl;
    file << "Convective Heat Transfer Coefficient :\t" << Convective_Heat_Transfer_Coefficient << "\t W / m2-K" << std::endl;

    file << std::endl;
}

long double Combustion_Pellet::Get_Pre_Heat_Zone_Transient_Term_Coeff()
{
    return Pre_Heat_Zone.Get_Density() * Pre_Heat_Zone.Get_Heat_Capacity();
}

long double Combustion_Pellet::Get_Reaction_Zone_Transient_Term_Coeff()
{
    return Reaction_Zone.Get_Density() * Reaction_Zone.Get_Heat_Capacity();
}

long double Combustion_Pellet::Get_Post_Combustion_Zone_Transient_Term_Coeff()
{
    return Post_Combustion_Zone.Get_Density() * Post_Combustion_Zone.Get_Heat_Capacity();
}

long double Combustion_Pellet::Get_Pre_Heat_Zone_Diffusion_Term_Coeff()
{
    return Pre_Heat_Zone.Get_Heat_Conductivity();
}

long double Combustion_Pellet::Get_Reaction_Zone_Diffusion_Term_Coeff()
{
    return Reaction_Zone.Get_Heat_Conductivity();
}

long double Combustion_Pellet::Get_Post_Combustion_Zone_Diffusion_Term_Coeff()
{
    return Post_Combustion_Zone.Get_Heat_Conductivity();
}

long double Combustion_Pellet::Get_Lateral_Surface_Heat_Loss_Term(long double T)
{
    return (4.0l / Diameter) * (
        Convective_Heat_Transfer_Coefficient * (T - Ambient_Temperature) + 
        Radiative_Emissivity * Stefan_Boltzmann_Constant * (std::pow(T, 4) - std::pow(Ambient_Temperature, 4))
    );
}

long double Combustion_Pellet::Get_Cross_Section_Area()
{
    return M_PI * Diameter * Diameter / 4.0l;
}

long double Combustion_Pellet::Get_Flat_Surface_Heat_Loss_Term(long double T)
{
    return Get_Cross_Section_Area() * (
        Convective_Heat_Transfer_Coefficient * (T - Ambient_Temperature) +
        Radiative_Emissivity * Stefan_Boltzmann_Constant * (std::pow(T, 4) - std::pow(Ambient_Temperature, 4))
    );
}