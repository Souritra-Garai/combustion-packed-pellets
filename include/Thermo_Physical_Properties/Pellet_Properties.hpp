#ifndef __PELLET_PROPS__
#define __PELLET_PROPS__

#include "Thermo_Physical_Properties/Coated_Particle.hpp"

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
            long double Packing_Volume_Fraction
        );

        long double Get_Density();
        long double Get_Heat_Capacity();
        long double Get_Heat_Conductivity();
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
            long double Packing_Volume_Fraction
        );

        long double Get_Density();
        long double Get_Heat_Capacity();
        long double Get_Heat_Conductivity();
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

long double Calc_Heat_Conductivity_Maxwell_Eucken_Bruggemann_Model(
    long double Particle_Volume_Fraction,
    long double Particle_Volume_Fraction_Maxwell_Eucken_Structure,
    long double Particle_Heat_Conductivity,
    long double Degassed_Fluid_Heat_Conductivity
);

#endif