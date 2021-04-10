#ifndef __COATED_PARTICLE__
#define __COATED_PARTICLE__

#include "Substance.hpp"

#include <math.h>

class Coated_Particle
{
    private:

        Substance *Core;
        Substance *Shell;

        long double Core_Diameter;
        long double Overall_Diameter;
        
    public:

        Coated_Particle(Substance CORE, Substance SHELL, long double CORE_DIAMETER, long double OVERALL_DIAMETER);

        long double Get_Density();
        long double Get_Heat_Capacity();
        long double Get_Heat_Conductivity();
};

#endif