#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <string.h>

#include "Combustion_Problem.hpp"

#define MAX_TIME_ITER   10

#define MOL_WT_Al   26.98154E-3 // kg / mol
#define MOL_WT_Ni   58.69E-3    // kg / mol
#define MOL_WT_Ar   39.9E-3     // kg / mol

// Declaring objects for holding
// Density (kg/m3), Heat Capacity (J/kg-K), Heat Conductivity (W/m-K), Molecular Weight (kg/mol)
// of pure substances at respective temperature zones
Substance Preheat_Zone_Aluminium(2700, 1060, 220, MOL_WT_Al);
Substance Preheat_Zone_Nickel(8908, 440, 66, MOL_WT_Ni);
Substance Preheat_Zone_Argon(0.89, 520, 0.016, MOL_WT_Ar);

Substance Post_Combustion_Zone_NiAl(5900, 717, 115, MOL_WT_Al + MOL_WT_Ni);
Substance Post_Combustion_Zone_Argon(0.37, 520, 0.055, MOL_WT_Ar);

Substance Reaction_Zone_Aluminium(2700, 1176, 130, MOL_WT_Al);
Substance Reaction_Zone_Nickel(8908, 670, 80, MOL_WT_Ni);
Substance Reaction_Zone_Argon(0.37, 520, 0.055, MOL_WT_Ar);
Substance Reaction_Zone_NiAl(5900, 717, 115, MOL_WT_Al + MOL_WT_Ni);

// Diameters of Ni Coated Al Particle
long double Reaction_Particle_Core_Diameter     = 65E-6; // m
long double Reaction_Particle_Overall_Diameter  = 79E-6; // m

// Objects to hold properties of Reaction Particles
// at respective temperature zones
Coated_Particle Pre_Heat_Zone_Ni_Coated_Al_Particle(
    Preheat_Zone_Aluminium,             // Core Substance
    Preheat_Zone_Nickel,                // Shell Substance
    Reaction_Particle_Core_Diameter,    // Core Diameter
    Reaction_Particle_Overall_Diameter  // Overall Diameter
);

Coated_Particle Reaction_Zone_Ni_Coated_Al_Particle(
    Reaction_Zone_Aluminium,            // Core Substance
    Reaction_Zone_Nickel,               // Shell Substance
    Reaction_Particle_Core_Diameter,    // Core Diameter
    Reaction_Particle_Overall_Diameter  // Shell Diameter
);

// Moles of Al and Ni in the Ni Coated Al Particle
long double Moles_Al_Reaction_Particle = Reaction_Zone_Ni_Coated_Al_Particle.Get_Moles_of_Core_Material();
long double Moles_Ni_Reaction_Particle = Reaction_Zone_Ni_Coated_Al_Particle.Get_Moles_of_Shell_Material();

// Limiting and Excess Agents
long double Moles_Limiting_Agent = std::min(Moles_Al_Reaction_Particle, Moles_Ni_Reaction_Particle);
long double Moles_Excess_Agent   = std::abs(Moles_Ni_Reaction_Particle - Moles_Al_Reaction_Particle);
long double Concentration_Limiting_Agent = std::min(
    Reaction_Zone_Ni_Coated_Al_Particle.Get_Core_Material_Concentration(),
    Reaction_Zone_Ni_Coated_Al_Particle.Get_Shell_Material_Concentration()
);
// Substance in excess in the Reaction Particle
Substance Excess_Agent = Moles_Al_Reaction_Particle < Moles_Ni_Reaction_Particle ? Reaction_Zone_Nickel : Reaction_Zone_Aluminium;

// Diameters of Product Particle
long double Product_Particle_Core_Diameter = Calc_Product_Particle_Core_Diameter(
    Moles_Excess_Agent, // Moles in excess of the excess agent
    Excess_Agent        // Excess agent substance
);

long double Product_Particle_Overall_Diameter = Calc_Product_Particle_Overall_Diameter(
    Moles_Limiting_Agent,           // Moles of Product produced
    Reaction_Zone_NiAl,             // Product substance
    Product_Particle_Core_Diameter  // Core diameter of Product Particle
);

// Objects to hold properties of Product Particles
// at respective temperature zones
Coated_Particle Reaction_Zone_NiAl_Particle(
    Excess_Agent,                       // Core Substance
    Reaction_Zone_NiAl,                 // Shell Substance
    Product_Particle_Core_Diameter,     // Core Diameter
    Product_Particle_Overall_Diameter   // Overall Diameter
);

Coated_Particle Post_Combustion_Zone_NiAl_Particle(
    Excess_Agent,                       // Core Substance
    Post_Combustion_Zone_NiAl,          // Shell Substance
    Product_Particle_Core_Diameter,     // Core Diameter
    Product_Particle_Overall_Diameter   // Overall Diameter
);

// Volume Fraction of Particles in the Pellet
long double Particle_Volume_Fraction = 0.5;

// Objects to hold Pellet Properties 
// at respective temperature zones
Pellet_Properties Pre_Heat_Zone_Pellet(
    Pre_Heat_Zone_Ni_Coated_Al_Particle,    // Particle
    Preheat_Zone_Argon,                     // Degassed Fluid Substance
    Particle_Volume_Fraction                // Volume Fraction Occupied by Particle
);

Pellet_Properties Post_Combustion_Zone_Pellet(
    Post_Combustion_Zone_NiAl_Particle,     // Particle
    Post_Combustion_Zone_Argon,             // Degassed Fluid Substance
    Particle_Volume_Fraction                // Volume Fraction Occupied by Particle
);

Reaction_Zone_Pellet_Properties Reaction_Zone_Pellet(
    Reaction_Zone_Ni_Coated_Al_Particle,    // Particle A
    Reaction_Zone_NiAl_Particle,            // Particle B
    Reaction_Zone_Argon,                    // Degassed Fluid Substance
    Particle_Volume_Fraction                // Volume Fraction Occupied by Particle
);

// Length and diameter of the combustion pellet
long double Pellet_Length = 6.35E-3;    // m
long double Pellet_Diameter = 6.35E-3;  // m

Kinetics Sundaram_et_al(
    "Sundaram et al 2013",
    465.23l * Concentration_Limiting_Agent,
    34.7E3l, 
    0
);

Reaction Combustion_Reaction(
    Sundaram_et_al,
    -118.4E3l,
    Concentration_Limiting_Agent,
    Particle_Volume_Fraction
);

// Matrix for Temperature at the grid points
// for each time point
std::vector<std::vector<long double>> T_MATRIX;

// Matrix for Conversion at the grid points
// for each time point
std::vector<std::vector<long double>> ETA_MATRIX;

unsigned int No_Grid_Points = 11;
long double Time_Step_Size = 0.2E-3;

int main(int argc, char** argv)
{   
    
    // Abstract object to represent all properties
    // of the combustion pellet
    Combustion_Pellet Pellet(
        Pellet_Diameter,                // Diameter of the Pellet
        Pellet_Length,                  // Length of the Pellet
        Pre_Heat_Zone_Pellet,           // Pellet Properties in Pre-Heat Zone
        Post_Combustion_Zone_Pellet,    // Pellet Properties in Post-Combustion Zone
        Reaction_Zone_Pellet            // Reaction Zone Pellet Properties
    );
    // Setting the environment of the Combustion Pellet
    Pellet.Set_Ambient_Temperature(298);    // K
    Pellet.Set_Convective_Heat_Transfer_Coefficient(19.68); // W/m2-K
    Pellet.Set_Radiative_Emissivity(0.25);  // Dimensionless

    Combustion_Problem my_Problem(
        No_Grid_Points,
        Time_Step_Size,
        Pellet,
        Combustion_Reaction,
        933,
        1911
    );

    std::cout << "Initialised Problem" << std::endl;

    std::vector<long double> Temperature_Array(No_Grid_Points, 298);

    std::vector<long double> Conversion_Array(No_Grid_Points, 0);

    Temperature_Array[0] = 1911;
    Conversion_Array[0] = 1;

    for (int i = 1; i < 1; i++)
    {
        Temperature_Array[i] = 933;
        Conversion_Array[i] = 0;
    }

    my_Problem.Set_Initial_Conditions(Temperature_Array, Conversion_Array);

    T_MATRIX.push_back(Temperature_Array);
    ETA_MATRIX.push_back(Conversion_Array);

    std::cout << "Solver set up." << std::endl;
    
    int n = 0;

    do
    {
        std::cout << "Time step : " << (n+1) * Time_Step_Size << " s" << std::endl;

        std::pair<std::vector<long double>, std::vector<long double>> Buffer = my_Problem.Iterate();

        T_MATRIX.push_back(Buffer.first);
        ETA_MATRIX.push_back(Buffer.second);

        n++;

        std::cout << "Completed " << n << " iterations\n" << std::endl;

    } while (n < MAX_TIME_ITER && my_Problem.Combustion_Not_Completed());

    std::cout << "Saving to file...\n";

    time_t rawtime;
    struct tm * timeinfo;
    char folder[100] = "solutions/";
    char T_file[100], Eta_file[100], config_file[100];

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );

    strcat(folder, asctime(timeinfo));
    strcpy(T_file, folder);
    strcpy(Eta_file, folder);
    strcpy(config_file, folder);

    mkdir(folder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    strcat(T_file, "/Temperature.csv");
    std::ofstream file(T_file);

    for (auto T_VECTOR = T_MATRIX.begin(); T_VECTOR < T_MATRIX.end(); T_VECTOR++)
    {
        for (auto T = T_VECTOR->begin(); T < T_VECTOR->end(); T++)
        {
            file << std::to_string((*T)) << ',';
        }

        file << std::endl;
    }

    file.close();

    strcat(config_file, "/Combustion_Config.txt");
    std::ofstream txt_file(config_file);

    txt_file << "Number of Time Steps :\t" << n+1 << std::endl;

    my_Problem.Write_to_File(txt_file);

    Pellet.Write_to_File(txt_file, "Ni Coated Al Pellet degassed with Ar");

    txt_file << "Pellet Volume Fraction Occupied by Solid Particles :\t" << Particle_Volume_Fraction << std::endl;

    txt_file << std::endl;

    Pre_Heat_Zone_Ni_Coated_Al_Particle.Write_to_File(txt_file, "Pre-Heat Zone Ni Coated Al Particle");
    Reaction_Zone_Ni_Coated_Al_Particle.Write_to_File(txt_file, "Reaction Zone Ni Coated Al Particle");
    Reaction_Zone_NiAl_Particle.Write_to_File(txt_file, "Reaction Zone Product NiAl Particle");
    Post_Combustion_Zone_NiAl_Particle.Write_to_File(txt_file, "Post-Combustion Zone Product NiAl Particle");

    Preheat_Zone_Aluminium.Write_to_File(txt_file, "Pre-Heat Zone Al");
    Preheat_Zone_Nickel.Write_to_File(txt_file, "Pre-Heat Zone Ni");
    Preheat_Zone_Argon.Write_to_File(txt_file, "Pre-Heat Zone Ar");

    Reaction_Zone_Aluminium.Write_to_File(txt_file, "Reaction Zone Al");
    Reaction_Zone_Nickel.Write_to_File(txt_file, "Reaction Zone Ni");
    Reaction_Zone_NiAl.Write_to_File(txt_file, "Reaction Zone NiAl");
    Reaction_Zone_Argon.Write_to_File(txt_file, "Reaction Zone Ar");

    Post_Combustion_Zone_NiAl.Write_to_File(txt_file, "Post-Combustion Zone NiAl");
    Post_Combustion_Zone_Argon.Write_to_File(txt_file, "Post-Combustion Zone Ar");

    txt_file.close();

    std::cout << "Solution saved to file.." << std::endl;

    return 0;    
}
