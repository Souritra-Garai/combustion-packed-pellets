#ifndef __PELLET_PROPS__
#define __PELLET_PROPS__

#include <fstream>

#include "Thermo_Physical_Properties/Coated_Particle.hpp"

// Volume fraction of particles in 
// Maxwell Eucken Structure
#define ALPHA_P_ME 0.7l

class Pellet_Properties
{
    private:
        
        long double Density;
        long double Heat_Capacity;
        long double Heat_Conductivity;

    public:

        Pellet_Properties(
            Coated_Particle Particle,
            Substance Degassed_Fluid,
            long double Packing_Volume_Fraction,
            long double (*Heat_Capacity_Calculating_Function) (
                long double Particle_Volume_Fraction,
                long double Particle_Heat_Conductivity,
                long double Degassed_Fluid_Heat_Conductivity
            )
        );

        long double Get_Density();
        long double Get_Heat_Capacity();
        long double Get_Heat_Conductivity();

        void Write_to_File(std::ofstream &File_Handler, const char *Name);
};

class Reaction_Zone_Pellet_Properties
{
    private:
        
        long double Density;
        long double Heat_Capacity;
        long double Heat_Conductivity;

    public:

        Reaction_Zone_Pellet_Properties(
            Coated_Particle Reactant,
            Coated_Particle Product,
            Substance Degassed_Fluid,
            long double Packing_Volume_Fraction,
            long double (*Heat_Capacity_Calculating_Function) (
                long double Particle_Volume_Fraction,
                long double Particle_Heat_Conductivity,
                long double Degassed_Fluid_Heat_Conductivity
            )
        );

        long double Get_Density();
        long double Get_Heat_Capacity();
        long double Get_Heat_Conductivity();

        void Write_to_File(std::ofstream &File_Handler, const char *Name);
};

long double Calc_Density(
    long double Particle_Volume_Fraction,
    long double Particle_Density,
    long double Degassed_Fluid_Density
);

long double Calc_Heat_Capacity(
    long double Particle_Volume_Fraction,
    long double Particle_Density,
    long double Particle_Heat_Capacity,
    long double Degassed_Fluid_Density,
    long double Degassed_Fluid_Heat_Capacity
);

long double Calc_Heat_Conductivity_Bruggemann_Model(
    long double Particle_Volume_Fraction,
    long double Particle_Heat_Conductivity,
    long double Degassed_Fluid_Heat_Conductivity
);

long double Calc_Heat_Conductivity_Maxwell_Eucken_Model(
    long double Particle_Volume_Fraction,
    long double Particle_Heat_Conductivity,
    long double Degassed_Fluid_Heat_Conductivity
);

long double Calc_Heat_Conductivity_Maxwell_Eucken_Bruggemann_Model(
    long double Particle_Volume_Fraction,
    long double Particle_Heat_Conductivity,
    long double Degassed_Fluid_Heat_Conductivity
);

#endif