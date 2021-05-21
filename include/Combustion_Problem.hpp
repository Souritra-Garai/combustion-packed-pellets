#ifndef __COMBUSTION_PROBLEM__
#define __COMBUSTION_PROBLEM__

#include <vector>
#include <fstream>

#include "Thermo_Physical_Properties/Combustion_Pellet.hpp"
#include "Kinetics/Reaction.hpp"

class Combustion_Problem
{
    private:

        unsigned int N;

        long double Delta_x;
        long double Delta_t;
        
        Combustion_Pellet Pellet;
        Reaction Combustion_Reaction;

        long double Ignition_Temperature;
        long double Adiabatic_Combustion_Temperature;

        long double *T_VECTOR;
        
        long double *ETA_VECTOR;

        bool In_Reaction_Zone(long double Temperature);
        bool In_Post_Combustion_Zone(long double Conversion, long double Temperature);

    public:

        Combustion_Problem(
            unsigned int No_of_Grid_Points,
            long double Time_Step_Size,
            Combustion_Pellet Pellet,
            Reaction Combustion_Reaction,
            long double Ignition_Temperature,
            long double Adiabatic_Combustion_Temperature
        );

        ~Combustion_Problem();

        void Set_Initial_Conditions(
            std::vector<long double> Temperature_Array,
            std::vector<long double> Conversion_Array
        );

        void Solve_and_Update_State();
        
        std::pair<std::vector<long double>, std::vector<long double>> Iterate();

        bool Combustion_Not_Completed();

        void Write_to_File(std::ofstream &File);
};

#endif