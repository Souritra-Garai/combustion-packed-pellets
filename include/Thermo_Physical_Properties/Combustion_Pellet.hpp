#ifndef __COMBUSTION_PELLET__
#define __COMBUSTION_PELLET__

#include <vector>
#include <fstream>

#include "Thermo_Physical_Properties/Pellet_Properties.hpp"

#define Stefan_Boltzmann_Constant 5.670374419E-8

class Combustion_Pellet
{
    private:
    
        long double Diameter;
        long double Length;

        Pellet_Properties Pre_Heat_Zone;
        Pellet_Properties Post_Combustion_Zone;
        Reaction_Zone_Pellet_Properties Reaction_Zone;

        long double Convective_Heat_Transfer_Coefficient;
        long double Radiative_Emissivity;
        long double Ambient_Temperature;

        long double Ignition_Temperature;

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
        void Set_Ignition_Temperature(long double Ignition_Temperature);

        long double Get_Pellet_Length();
        long double Get_Cross_Section_Area();
        long double Get_Ignition_Temperature();

        long double Get_Pre_Heat_Zone_Transient_Term_Coeff();
        long double Get_Reaction_Zone_Transient_Term_Coeff();
        long double Get_Post_Combustion_Zone_Transient_Term_Coeff();

        long double Get_Pre_Heat_Zone_Diffusion_Term_Coeff();
        long double Get_Reaction_Zone_Diffusion_Term_Coeff();
        long double Get_Post_Combustion_Zone_Diffusion_Term_Coeff();

        long double Get_Lateral_Surface_Heat_Loss_Term(long double Temperature);
        long double Get_Flat_Surface_Heat_Loss_Term(long double Temperature);

        void Write_to_File(std::ofstream &File_Handler, const char *Name);
};

#endif