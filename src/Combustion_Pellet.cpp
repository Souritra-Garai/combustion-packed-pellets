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

std::pair<long double, long double> Combustion_Pellet::Get_Heat_Loss_Polynomial(long double T)
{
    long double conv_coeff = 4 * Convective_Heat_Transfer_Coefficient / Diameter;
    long double rad_coeff  = 4 * Radiative_Emissivity * Stefan_Boltzmann_Constant / Diameter;

    return std::pair<long double, long double> (
        conv_coeff + 4 * rad_coeff * pow(T, 3),
        conv_coeff * Ambient_Temperature + rad_coeff * (3 * pow(T, 4) + pow(Ambient_Temperature, 4))
    );
}

void Combustion_Pellet::Setup_Pre_Heat_Zone_Equation(
    long double &E,
    long double &F,
    long double &G,
    long double &B,
    long double T,
    long double Dx,
    long double Dt
)
{
    E = G = - Pre_Heat_Zone.Get_Heat_Conductivity() / (Dx * Dx);

    std::pair<long double, long double> polynomial = Get_Heat_Loss_Polynomial(T);

    F = Pre_Heat_Zone.Get_Density() * Pre_Heat_Zone.Get_Heat_Capacity() / Dt + 2 * Pre_Heat_Zone.Get_Heat_Conductivity() / (Dx * Dx) + polynomial.first;
    B = Pre_Heat_Zone.Get_Density() * Pre_Heat_Zone.Get_Heat_Capacity() * T / Dt + polynomial.second;
}

void Combustion_Pellet::Setup_Reaction_Zone_Equation(
    long double &E,
    long double &F,
    long double &G,
    long double &B,
    long double T,
    long double Dx,
    long double Dt
)
{
    E = G = - Reaction_Zone.Get_Heat_Conductivity() / (Dx * Dx);

    std::pair<long double, long double> polynomial = Get_Heat_Loss_Polynomial(T);

    F = Reaction_Zone.Get_Density() * Reaction_Zone.Get_Heat_Capacity() / Dt + 2 * Reaction_Zone.Get_Heat_Conductivity() / (Dx * Dx) + polynomial.first;
    B = Reaction_Zone.Get_Density() * Reaction_Zone.Get_Heat_Capacity() * T / Dt + polynomial.second;
}

void Combustion_Pellet::Setup_Post_Combustion_Zone_Equation(
    long double &E,
    long double &F,
    long double &G,
    long double &B,
    long double T,
    long double Dx,
    long double Dt
)
{
    E = G = - Post_Combustion_Zone.Get_Heat_Conductivity() / (Dx * Dx);

    std::pair<long double, long double> polynomial = Get_Heat_Loss_Polynomial(T);

    F = Post_Combustion_Zone.Get_Density() * Post_Combustion_Zone.Get_Heat_Capacity() / Dt + 2 * Post_Combustion_Zone.Get_Heat_Conductivity() / (Dx * Dx) + polynomial.first;
    B = Post_Combustion_Zone.Get_Density() * Post_Combustion_Zone.Get_Heat_Capacity() * T / Dt + polynomial.second;
}

void Combustion_Pellet::Setup_X0_Isothermal_BC_Equation(
    long double &E,
    long double &F,
    long double &G,
    long double &B,
    long double T,
    long double Dx,
    long double Dt
)
{
    E = 0;
    F = 1;
    G = 0;
    B = T;
}


void Combustion_Pellet::Setup_X0_Adiabatic_Wall_BC_Equation(
    long double &E,
    long double &F,
    long double &G,
    long double &B,
    long double T,
    long double Dx,
    long double Dt
)
{
    E = 0;
    F = 1;
    G = -1;
    B = 0;
}

void Combustion_Pellet::Setup_X0_Ambient_Heat_Loss_BC_Equation(
    long double &E,
    long double &F,
    long double &G,
    long double &B,
    long double T,
    long double Dx,
    long double Dt,
    int z
)
{
    long double lambda = z * (z-1) * Post_Combustion_Zone.Get_Heat_Conductivity() / 2 + z * (z-2) * Reaction_Zone.Get_Heat_Conductivity() / (-1) + (z-1) * (z-2) * Pre_Heat_Zone.Get_Heat_Conductivity() / 2; 
    E = 0;
    F = - (lambda / Dx) - Convective_Heat_Transfer_Coefficient - 4 * Radiative_Emissivity * Stefan_Boltzmann_Constant * pow(T, 3);
    G = lambda / Dx;
    B = - Convective_Heat_Transfer_Coefficient * Ambient_Temperature - Radiative_Emissivity * Stefan_Boltzmann_Constant * (3*pow(T, 4) + pow(Ambient_Temperature, 4));
}

void Combustion_Pellet::Setup_XM_Adiabatic_Wall_BC_Equation(
    long double &E,
    long double &F,
    long double &G,
    long double &B,
    long double T,
    long double Dx,
    long double Dt
)
{
    E = -1;
    F = 1;
    G = 0;
    B = 0;
}

void Combustion_Pellet::Setup_XM_Ambient_Heat_Loss_BC_Equation(
    long double &E,
    long double &F,
    long double &G,
    long double &B,
    long double T,
    long double Dx,
    long double Dt,
    int z
)
{
    long double lambda = z * (z-1) * Post_Combustion_Zone.Get_Heat_Conductivity() / 2 + z * (z-2) * Reaction_Zone.Get_Heat_Conductivity() / (-1) + (z-1) * (z-2) * Pre_Heat_Zone.Get_Heat_Conductivity() / 2; 
    
    E = lambda / Dx;
    F = - (lambda / Dx) - Convective_Heat_Transfer_Coefficient - 4 * Radiative_Emissivity * Stefan_Boltzmann_Constant * pow(T, 3);
    G = 0;
    B = - Convective_Heat_Transfer_Coefficient * Ambient_Temperature - Radiative_Emissivity * Stefan_Boltzmann_Constant * (3*pow(T, 4) + pow(Ambient_Temperature, 4));
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

void Write_to_File(std::ofstream &file, const char *name, Substance s)
{
    file << std::endl << name << " Properties" << std::endl;

    file << "Density :\t" << s.Get_Density() << "\t kg / m3" << std::endl;
    file << "Heat Capacity :\t" << s.Get_Heat_Capacity() << "\t J / kg - K" << std::endl;
    file << "Heat Conductivity :\t" << s.Get_Heat_Conductivity() << "\t W / m - K" << std::endl;

    file << std::endl;
}

void Write_to_File(std::ofstream &file, const char *name, Coated_Particle s)
{
    file << std::endl << name << " Properties" << std::endl;

    file << "Density :\t" << s.Get_Density() << "\t kg / m3" << std::endl;
    file << "Heat Capacity :\t" << s.Get_Heat_Capacity() << "\t J / kg - K" << std::endl;
    file << "Heat Conductivity :\t" << s.Get_Heat_Conductivity() << "\t W / m - K" << std::endl;

    file << std::endl;
}

