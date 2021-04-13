#include <iostream>
#include <fstream>
#include <vector>

#include "Combustion_Problem.hpp"

#define MAX_TIME_ITER   1000
#define NO_GRID_POINTS  10

#define TIME_STEP_SIZE  1E-4

#define phi_P   0.5

#define MOL_WT_Al   26.98154E-3 // kg / mol
#define MOL_WT_Ni   58.69E-3    // kg / mol

Substance Preheat_Zone_Aluminium(2700, 1060, 220);
Substance Preheat_Zone_Nickel(8908, 440, 66);
Substance Preheat_Zone_Argon(0.89, 520, 0.016);

Substance Post_Combustion_Zone_NiAl(5900, 717, 115);
Substance Post_Combustion_Zone_Argon(0.37, 520, 0.055);

Substance Reaction_Zone_Aluminium(2700, 1176, 130);
Substance Reaction_Zone_Nickel(8908, 670, 80);
Substance Reaction_Zone_Argon(0.37, 520, 0.055);
Substance Reaction_Zone_NiAl(5900, 717, 115);

long double Reaction_Particle_Core_Diameter     = 65E-6; // m
long double Reaction_Particle_Overall_Diameter  = 79E-6; // m

long double Mol_Ni = (M_PIl / 6) * (
    pow(Reaction_Particle_Overall_Diameter, 3) -
    pow(Reaction_Particle_Core_Diameter, 3)
) * Preheat_Zone_Nickel.Get_Density() / MOL_WT_Ni;

long double Mol_Al = (M_PIl / 6) *
    pow(Reaction_Particle_Core_Diameter, 3) *
    Preheat_Zone_Aluminium.Get_Density() /
    MOL_WT_Al;

long double Product_Particle_Core_Diameter = pow(
    (Mol_Ni - Mol_Al) * MOL_WT_Ni /
    (Reaction_Zone_Nickel.Get_Density() * (M_PIl / 6)),
    1.0/3.0
);

long double Product_Particle_Overall_Diameter = pow(
    pow(Product_Particle_Core_Diameter, 3) +
    Mol_Al * (MOL_WT_Al + MOL_WT_Ni) /
    (Reaction_Zone_NiAl.Get_Density() * (M_PIl / 6)),
    1.0/3.0
);

Coated_Particle Pre_Heat_Zone_Ni_Coated_Al_Particle(
    Preheat_Zone_Aluminium,
    Preheat_Zone_Nickel,
    Reaction_Particle_Core_Diameter,
    Reaction_Particle_Overall_Diameter
);

Coated_Particle Reaction_Zone_Ni_Coated_Al_Particle(
    Reaction_Zone_Aluminium,
    Reaction_Zone_Nickel,
    Reaction_Particle_Core_Diameter,
    Reaction_Particle_Overall_Diameter
);

Coated_Particle Reaction_Zone_NiAl_Particle(
    Reaction_Zone_Aluminium,
    Reaction_Zone_NiAl,
    Product_Particle_Core_Diameter,
    Product_Particle_Overall_Diameter
);

Coated_Particle Post_Combustion_Zone_Particle(
    Reaction_Zone_Aluminium,
    Post_Combustion_Zone_NiAl,
    Product_Particle_Core_Diameter,
    Product_Particle_Overall_Diameter
);

Pellet_Properties Pre_Heat_Zone_Pellet(
    Pre_Heat_Zone_Ni_Coated_Al_Particle,
    Preheat_Zone_Argon,
    phi_P
);

Pellet_Properties Post_Combustion_Zone_Pellet(
    Post_Combustion_Zone_Particle,
    Post_Combustion_Zone_Argon,
    phi_P
);

Reaction_Zone_Pellet_Properties Reaction_Zone_Pellet(
    Reaction_Zone_Ni_Coated_Al_Particle,
    Reaction_Zone_NiAl_Particle,
    Reaction_Zone_Argon,
    phi_P
);

Combustion_Pellet my_Pellet(
    6.35E-3,
    6.35E-3,
    Pre_Heat_Zone_Pellet,
    Post_Combustion_Zone_Pellet,
    Reaction_Zone_Pellet
);

Kinetics Sundaram_et_al(
    "Sundaram et al 2013",
    465l,
    34.7E3l, 
    0
);

Reaction Combustion_Reaction(
    Sundaram_et_al,
    - 50E3,
    Mol_Al / ((M_PIl / 6) * pow(Reaction_Particle_Overall_Diameter, 3)),
    phi_P
);

// Matrix for Temperature at the grid points
// for each time point
std::vector<std::vector<long double>> T_MATRIX;

// Matrix for Conversion at the grid points
// for each time point
std::vector<std::vector<long double>> ETA_MATRIX;

// Function to check Temperature vectors for two consecutive time steps
// are equal, that is steady state has been reached
bool has_converged(std::vector<long double>, std::vector<long double>);

int main(int argc, char** argv)
{   
    my_Pellet.Set_Ambient_Temperature(298);
    my_Pellet.Set_Convective_Heat_Transfer_Coefficient(19.68);
    my_Pellet.Set_Radiative_Emissivity(0.25);

    Combustion_Problem my_Problem(
        NO_GRID_POINTS,
        TIME_STEP_SIZE,
        my_Pellet,
        Combustion_Reaction,
        933,
        1911
    );

    std::cout << "Initialised Problem" << std::endl;

    std::vector<long double> Temperature_Array(NO_GRID_POINTS, 298);
    Temperature_Array[0] = 1911;

    for (int i = 0; i < 10; i++)

        Temperature_Array[i] = 933;

    std::vector<long double> Conversion_Array(NO_GRID_POINTS, 0);
    Conversion_Array[0] = 1;

    my_Problem.Set_Initial_Conditions(Temperature_Array, Conversion_Array);

    T_MATRIX.push_back(Temperature_Array);
    ETA_MATRIX.push_back(Conversion_Array);

    std::cout << "Solver set up." << std::endl;
    
    int n = 0;

    do
    {
        std::cout << "Time step : " << (n+1) * TIME_STEP_SIZE << " s" << std::endl;

        std::pair<std::vector<long double>, std::vector<long double>> Buffer = my_Problem.Iterate();

        T_MATRIX.push_back(Buffer.first);
        ETA_MATRIX.push_back(Buffer.second);

        n++;

        std::cout << "Completed " << n << " iterations\n" << std::endl;

    } while (n < MAX_TIME_ITER && !has_converged(T_MATRIX[n], T_MATRIX[n-1]));

    std::cout << "Saving to file...\n";

    std::ofstream file("solutions/solution_1.csv");

    for (auto T_VECTOR = T_MATRIX.begin(); T_VECTOR < T_MATRIX.end(); T_VECTOR++)
    {
        for (auto T = T_VECTOR->begin(); T < T_VECTOR->end(); T++)
        {
            file << std::to_string((*T)) << ',';
        }

        file << std::endl;
    }

    file.close();
    std::cout << "Solution saved to file.." << std::endl;

    return 0;    
}

bool is_close(long double a, long double b)
{
    if (-0.01 < (a - b) && (a - b) < 0.01) return true;

    else return false;
}

bool has_converged(std::vector<long double> VECTOR_X, std::vector<long double> VECTOR_Y)
{
    bool flag = true;
    std::vector<long double>::iterator X, Y;
    
    for (X = VECTOR_X.begin(), Y = VECTOR_Y.begin(); X < VECTOR_X.end() && flag; X++, Y++)

        flag = is_close(*X, *Y);

    return flag;
}

