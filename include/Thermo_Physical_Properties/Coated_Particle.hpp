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
        
    public:

        Coated_Particle(Substance CORE, Substance SHELL, long double CORE_DIAMETER, long double OVERALL_DIAMETER);

        long double Get_Density();
        long double Get_Heat_Capacity();
        long double Get_Heat_Conductivity();
};

#endif