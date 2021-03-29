#include <vector>
#include <math.h>

#include "TDMA_solver.hpp"

#define Stefan_Boltzmann_Constant 5.670374419E-8

class Temperature_Iterator
{
    private:

        unsigned int N;

        long double L;
        long double D;

        long double Delta_t;
        long double Delta_x;

        long double T_atm;
        long double T_final;

        long double* rho;
        long double* C_p;
        long double* lambda;

        long double h_conv;
        long double epsilon;

        std::vector<long double> T_VECTOR;
        std::vector<long double>::iterator T;

        std::vector<long double> E_VECTOR;
        std::vector<long double> F_VECTOR;
        std::vector<long double> G_VECTOR;
        std::vector<long double> B_VECTOR;

        std::vector<long double>::iterator E, F, G, B;

        TDMA_solver::solver SOLVER;

        long double Calc_E(unsigned int);
        long double Calc_G(unsigned int);
        long double Calc_F(unsigned int, long double);
        long double Calc_B(unsigned int, long double);

        void Reset_Banded_Matrix_Iterators();
        bool In_Range_Banded_Matrix_Iterators();
        void Increment_Banded_Matrix_Iterators();

        void Apply_BC_Banded_Matrix();
        void Apply_BC_B_Vector();

    public:
        
        Temperature_Iterator(unsigned int N);

        void Set_Temperatures               (long double Atmospheric_Temperature, long double Final_Temperature);
        void Set_Time_Step_Length           (long double Dt);
        void Set_Pellet_Dimensions          (long double Length,    long double Diameter);
        void Set_Curved_Surface_Heat_Losses (long double Convective_Heat_Transfer_Coefficient,  long double Emissivity);
        void Set_Thermophysical_Properties  (long double* Density,  long double* Heat_Capacity, long double* Heat_Diffusivity);

        void Apply_Initial_Condition(std::vector<long double> Initial_Temperature_Vector);
        
        void Setup_Matrix_Equation();
        
        std::vector<long double> Get_Solution();
};