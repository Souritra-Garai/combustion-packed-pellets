#ifndef __COMBUSTION_PROBLEM__
#define __COMBUSTION_PROBLEM__

#include <vector>
#include <fstream>

#include "Thermo_Physical_Properties/Combustion_Pellet.hpp"
#include "Kinetics/Reaction.hpp"
#include "TDMA_solver.hpp"

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

        std::vector<long double> T_VECTOR;
        std::vector<long double>::iterator T;

        std::vector<long double> ETA_VECTOR;
        std::vector<long double>::iterator ETA;

        std::vector<long double> E_VECTOR;
        std::vector<long double> F_VECTOR;
        std::vector<long double> G_VECTOR;
        std::vector<long double> B_VECTOR;

        std::vector<long double>::iterator E, F, G, B;

        TDMA_solver::solver SOLVER;

        void Reset_Equation_Iterators();
        bool In_Range_Equation_Iterators();
        void Increment_Equation_Iterators();

        bool In_Reaction_Zone(long double Temperature);
        bool In_Post_Combustion_Zone(long double Conversion, long double Temperature);

        void Setup_Solver();
        void Solve_and_Update_State();

    public:

        Combustion_Problem(
            unsigned int No_of_Grid_Points,
            long double Time_Step_Size,
            Combustion_Pellet Pellet,
            Reaction Combustion_Reaction,
            long double Ignition_Temperature,
            long double Adiabatic_Combustion_Temperature
        );

        void Set_Initial_Conditions(
            std::vector<long double> Temperature_Array,
            std::vector<long double> Conversion_Array
        );

        std::pair<std::vector<long double>, std::vector<long double>> Iterate();

        bool Combustion_Not_Completed();

        void Write_to_File(std::ofstream &File);
};

#endif