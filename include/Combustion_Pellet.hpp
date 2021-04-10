#ifndef __COMBUSTION_PELLET__
#define __COMBUSTION_PELLET__

#include <vector>

#include "Pellet_Properties.hpp"

#define Stefan_Boltzmann_Constant 5.670374419E-8

class Combustion_Pellet
{
    private:
    
        long double Diameter;
        long double Length;

        Pellet_Properties *Pre_Heat_Zone;
        Pellet_Properties *Post_Combustion_Zone;
        Reaction_Zone_Pellet_Properties *Reaction_Zone;

        long double Convective_Heat_Transfer_Coefficient;
        long double Radiative_Emissivity;
        long double Ambient_Temperature;

        std::pair<long double, long double> Get_Heat_Loss_Polynomial(long double Temperature);

    public:

        Combustion_Pellet(
            long double Diameter,
            long double Length,
            Pellet_Properties Pre_Heat_Zone,
            Pellet_Properties Post_Combustion_Zone,
            Reaction_Zone_Pellet_Properties Reaction_Zone
        );

        void Set_Convective_Heat_Transfer_Coefficient(long double Convective_Heat_Transfer_Coefficient);
        void Set_Radiative_Emissivity(long double Radiative_Emissivity);
        void Set_Ambient_Temperature(long double Ambient_Temperature);

        void Setup_Pre_Heat_Zone_Equation(
            long double &E,
            long double &F,
            long double &G,
            long double &B,
            long double T,
            long double Delta_x,
            long double Delta_t
        );
        void Setup_Reaction_Zone_Equation(
            long double &E,
            long double &F,
            long double &G,
            long double &B,
            long double T,
            long double Delta_x,
            long double Delta_t
        );
        void Setup_Post_Combustion_Zone_Equation(
            long double &E,
            long double &F,
            long double &G,
            long double &B,
            long double T,
            long double Delta_x,
            long double Delta_t
        );

        void Setup_X0_Isothermal_BC_Equation(
            long double &E,
            long double &F,
            long double &G,
            long double &B,
            long double T,
            long double Delta_x,
            long double Delta_t
        );
        void Setup_X0_Adiabatic_Wall_BC_Equation(
            long double &E,
            long double &F,
            long double &G,
            long double &B,
            long double T,
            long double Delta_x,
            long double Delta_t
        );
        void Setup_X0_Ambient_Heat_Loss_BC_Equation(
            long double &E,
            long double &F,
            long double &G,
            long double &B,
            long double T,
            long double Delta_x,
            long double Delta_t,
            int zone = 2
        );

        void Setup_XM_Adiabatic_Wall_BC_Equation(
            long double &E,
            long double &F,
            long double &G,
            long double &B,
            long double T,
            long double Delta_x,
            long double Delta_t
        );
        void Setup_XM_Ambient_Heat_Loss_BC_Equation(
            long double &E,
            long double &F,
            long double &G,
            long double &B,
            long double T,
            long double Delta_x,
            long double Delta_t,
            int zone = 2
        );
};

#endif