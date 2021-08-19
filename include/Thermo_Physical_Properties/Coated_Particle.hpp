#ifndef __COATED_PARTICLE__
#define __COATED_PARTICLE__

#include "Thermo_Physical_Properties/Substance.hpp"

#include <math.h>

class Coated_Particle
{
    private:

        long double Density;
        long double Heat_Capacity;
        long double Heat_Conductivity;

        long double Core_Diameter;
        long double Overall_Diameter;

        long double Moles_of_Core_Material;
        long double Moles_of_Shell_Material;
        
    public:

        Coated_Particle(Substance CORE, Substance SHELL, long double CORE_DIAMETER, long double OVERALL_DIAMETER);

        long double Get_Density();
        long double Get_Heat_Capacity();
        long double Get_Heat_Conductivity();

        void Write_to_File(std::ofstream &File_Handler, const char *Name);

        long double Get_Moles_of_Core_Material();
        long double Get_Moles_of_Shell_Material();
        long double Get_Core_Material_Concentration();
        long double Get_Shell_Material_Concentration();
};

long double Calc_Product_Particle_Core_Diameter(
    long double Moles_of_Core_Substance,
    Substance Core_Substance
);

long double Calc_Product_Particle_Overall_Diameter(
    long double Moles_of_Shell_Substance,
    Substance Shell_Substance,
    long double Core_Diameter
);

#endif