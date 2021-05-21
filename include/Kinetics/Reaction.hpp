#ifndef __REACTION__
#define __REACTION__

#include "Kinetics/Kinetics.hpp"
#include <vector>

class Reaction
{
    private :

        Kinetics Reaction_Kinetics;

        long double Enthalpy_of_Reaction;

        long double Initial_Concentration;

        long double Particle_Volume_Fraction;

    public :

        Reaction(
            Kinetics Reaction_Kinetics,
            long double Enthalpy_of_Reaction,
            long double Initial_Concentration,
            long double Particle_Volume_Fraction
        );

        std::pair<long double, long double> Get_Linear_Expression(
            long double Temperature,
            long double Conversion,
            long double Delta_t
        );

        long double Get_Omega(
            long double Temperature,
            long double Conversion
        );

        long double Get_Energy_Gen_Rate(
            long double Omega
        );

        void Update_Conversion(
            long double &Conversion,
            long double Temperature,
            long double Delta_t
        );
};

#endif