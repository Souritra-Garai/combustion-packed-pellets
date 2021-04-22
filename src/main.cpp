#include <iostream>
#include <vector>

#include "Combustion_Problem.hpp"
#include "File_Utilities.hpp"

#define MAX_TIME_ITER   10000
#define NO_GRID_POINTS  1001

#define IGNITION_LENGTH 0.5E-3  //  m

#define MOL_WT_Al   26.98154E-3 // kg / mol
#define MOL_WT_Ni   58.69E-3    // kg / mol
#define MOL_WT_Ar   39.9E-3     // kg / mol

#define ENTHALPY_CHANGE_COMBUSTION_REACTION  -118.4E3l   // J / mol.

#define HEAT_CONDUCTIVITY_CALC_FUNC Calc_Heat_Conductivity_Bruggemann_Model

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

// Length and diameter of the combustion pellet
long double Pellet_Length = 6.35E-3;    // m
long double Pellet_Diameter = 6.35E-3;  // m

// Object to model kinetics for the Ni-Al reaction
// Kinetics Sundaram_et_al(
//     "Sundaram et al 2013",  //  Name of Kinetics Model
//     465.23l,                //  Pre - Exponential Factor
//     34.7E3l,                //  Activation Energy
//     0,                      //  Reaction Order in (1 - eta)
//     0                       //  Reaction Order in eta
// );

// Object to model kinetics for the Ni-Al reaction
Kinetics Du_et_al(
    "Du et al 2013",  //  Name of Kinetics Model
    541.92l,                //  Pre - Exponential Factor
    26.0E3l,                //  Activation Energy
    0,                      //  Reaction Order in (1 - eta)
    0                       //  Reaction Order in eta
);

// Kinetics White_et_al(
//     "White et al 2009",     //  Name of Kinetics Model
//     3.7E19l,                //  Pre - Exponential Factor
//     347.29E3l,              //  Activation Energy
//     1,                      //  Reaction Order in (1 - eta)
//     0                       //  Reaction Order in eta
// );

// Kinetics Maiti_et_al(
//     "Maiti et al 2015",     //  Name of Kinetics Model
//     2.215E25,               //  Pre - Exponential Factor
//     448.4E3l,               //  Activation Energy
//     0.439,                  //  Reaction Order in (1 - eta)
//     1.5                     //  Reaction Order in eta
// );

// Matrix for Temperature at the grid points
// for each time point
std::vector<std::vector<long double>> T_MATRIX;

// Matrix for Conversion at the grid points
// for each time point
std::vector<std::vector<long double>> ETA_MATRIX;

unsigned int No_Grid_Points = NO_GRID_POINTS;
long double Time_Step_Size = 0.00002;

int main(int argc, char** argv)
{   
    // Volume Fraction of Particles in the Pellet
    long double Particle_Volume_Fraction = 0.5;

    // Objects to hold Pellet Properties 
    // at respective temperature zones
    Pellet_Properties Pre_Heat_Zone_Pellet(
        Pre_Heat_Zone_Ni_Coated_Al_Particle,    // Particle
        Preheat_Zone_Argon,                     // Degassed Fluid Substance
        Particle_Volume_Fraction,               // Volume Fraction Occupied by Particle
        HEAT_CONDUCTIVITY_CALC_FUNC             // Heat Capacity Calculating Function
    );

    Pellet_Properties Post_Combustion_Zone_Pellet(
        Post_Combustion_Zone_NiAl_Particle,     // Particle
        Post_Combustion_Zone_Argon,             // Degassed Fluid Substance
        Particle_Volume_Fraction,               // Volume Fraction Occupied by Particle
        HEAT_CONDUCTIVITY_CALC_FUNC             // Heat Capacity Calculating Function
    );

    Reaction_Zone_Pellet_Properties Reaction_Zone_Pellet(
        Reaction_Zone_Ni_Coated_Al_Particle,    // Particle A
        Reaction_Zone_NiAl_Particle,            // Particle B
        Reaction_Zone_Argon,                    // Degassed Fluid Substance
        Particle_Volume_Fraction,               // Volume Fraction Occupied by Particle
        HEAT_CONDUCTIVITY_CALC_FUNC             // Heat Capacity Calculating Function
    );
    
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

    // Interface between kinetics and combustion problem solver
    Reaction Combustion_Reaction(
        // Sundaram_et_al,                         // Kinetics model
        Du_et_al,                               // Kinetics model
        // White_et_al,                            // Kinetics model
        // Maiti_et_al,                            // Kinetics model
        ENTHALPY_CHANGE_COMBUSTION_REACTION,    // Enthalpy change for the reaction
        Concentration_Limiting_Agent,           // Concentration of the limiting agent
        Particle_Volume_Fraction                // Volume fraction occupied by Particle
    );

    Combustion_Problem Ni_Al_Pellet_Combustion(
        No_Grid_Points,         // Number of grid points in x-axis (including both boundary points)
        Time_Step_Size,         // Length of Time Steps
        Pellet,                 // Combustion_Pellet object holding pellet properties at different temperature zones
        Combustion_Reaction,    // Reaction object serving as interface between Kinetics and Combustion_Problem
        933,                    // Ignition temperature defining the reaction zone
        1911                    // Adiabatic Combustion Temperature
    );

    std::cout << "Initialised Problem" << std::endl;

    // Setting up the initial conditions
    std::vector<long double> Temperature_Array  (No_Grid_Points, 298);
    std::vector<long double> Conversion_Array   (No_Grid_Points, 0.000001);

    Temperature_Array[0] = 1911;
    Conversion_Array[0] = 0.999999;

    for (int i = 1; i * (Pellet_Length / (No_Grid_Points-1)) < IGNITION_LENGTH; i++)
    {
        Temperature_Array[i] = 1200;
        // Conversion_Array[i] = 0;
    }

    Ni_Al_Pellet_Combustion.Set_Initial_Conditions(Temperature_Array, Conversion_Array);

    T_MATRIX.push_back(Temperature_Array);
    ETA_MATRIX.push_back(Conversion_Array);

    std::cout << "Solver set up." << std::endl;
    
    int n = 0;

    do
    {
        std::cout << "Time step : " << (n+1) * Time_Step_Size << " s" << std::endl;

        std::pair<std::vector<long double>, std::vector<long double>> Buffer = Ni_Al_Pellet_Combustion.Iterate();

        T_MATRIX.push_back(Buffer.first);
        ETA_MATRIX.push_back(Buffer.second);

        n++;

        std::cout << "Completed " << n << " iterations\n" << std::endl;

    } while (n < MAX_TIME_ITER && Ni_Al_Pellet_Combustion.Combustion_Not_Completed());

    std::cout << "Saving to file...\n";

    Make_Solution_Folder();

    Save_Temperature_Data(T_MATRIX);
    Save_Conversion_Data(ETA_MATRIX);
    Save_Combustion_Config(n, Particle_Volume_Fraction, Pellet, Ni_Al_Pellet_Combustion);    

    std::cout << "Solution saved to file.." << std::endl;

    return 0;    
}
